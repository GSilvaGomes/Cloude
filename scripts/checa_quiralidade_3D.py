"""Confere a estereoquimica de estruturas 3D contra os SMILES de referencia.

Serve para validar o que saiu de um gerador de 3D (por exemplo o tradutor web
que o ChimeraX usa no `open smiles:`) antes de usar as estruturas em docking:
se um centro quiral vier invertido, o residuo deixa de ser L e a molecula e
outra, com pose diferente no sitio. Minimizar depois nao corrige isso --
a minimizacao preserva a quiralidade que recebeu.

COMO A CHECAGEM E FEITA

A quiralidade sai da geometria, nao das ordens de ligacao: para cada centro,
calcula-se o volume com sinal dos tres vizinhos pesados em torno do atomo
central, e compara-se o sinal com o de uma referencia 3D construida a partir
do SMILES (que e L por construcao). Sinais iguais, mesma configuracao.

Esse caminho foi escolhido porque PDB nao guarda ordem de ligacao. Tentar
reconstrui-la falha nestas moleculas: o RDKit marca os aneis de cinco
membros (pGlu, prolina) como aromaticos e nao consegue kekulizar. Atribuir
CIP a partir de ordens de ligacao erradas daria um resultado sem valor, ja
que a prioridade CIP depende delas. O volume com sinal nao depende.

Uso:
    python3 scripts/checa_quiralidade_3D.py arquivo.pdb ...
    python3 scripts/checa_quiralidade_3D.py estruturas_3D/

O nome do arquivo e casado com o registro em Smiles_mcro_analogos.txt.
Nomes como `L-mcro3` casam com `mcr3`; use --ref NOME para forcar.
"""
import os
import re
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog('rdApp.*')

RAIZ = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
REF = os.path.join(RAIZ, "Smiles_mcro_analogos.txt")


def carrega_referencia(path):
    ref, nome, linhas = {}, None, []
    for raw in open(path):
        s = raw.strip()
        if not s or s.startswith("#"):
            continue
        if s.startswith(">"):
            nome = s[1:].split("|")[0].strip()
            linhas = []
            ref[nome] = linhas
        elif nome:
            linhas.append(s)
    return {k: {"codigo": v[1], "smiles": v[2]}
            for k, v in ref.items() if len(v) >= 3}


def nome_para_ref(base, ref):
    """L-mcro3 -> mcr3 ; L-mcro-Y / L-mcro -> mcro."""
    if base in ref:
        return base
    b = base.lower().replace("l-", "").replace("_", "-")
    b = re.sub(r"-(y|min|em|opt)$", "", b)
    m = re.match(r"^mcro?(\d+)$", b)
    if m and f"mcr{m.group(1)}" in ref:
        return f"mcr{m.group(1)}"
    if b in ("mcro", "mcr0") and "mcro" in ref:
        return "mcro"
    return b if b in ref else None


def le_3d(path):
    """Le sem sanitizar: so precisamos de elementos, conectividade e xyz."""
    ext = os.path.splitext(path)[1].lower()
    if ext in (".pdb", ".ent"):
        m = Chem.MolFromPDBFile(path, removeHs=False, sanitize=False)
    elif ext in (".sdf", ".mol"):
        m = next(iter(Chem.SDMolSupplier(path, removeHs=False, sanitize=False)), None)
    elif ext == ".mol2":
        m = Chem.MolFromMol2File(path, removeHs=False, sanitize=False)
    else:
        return None, f"extensao nao suportada ({ext})"
    if m is None:
        return None, "nao foi possivel ler o arquivo"
    if m.GetNumConformers() == 0:
        return None, "sem coordenadas 3D"
    return m, None


def sinal_volume(conf, centro, vizinhos):
    """Sinal do produto misto dos vetores centro->vizinho. Define a mao."""
    c = conf.GetAtomPosition(centro)
    v = [conf.GetAtomPosition(i) - c for i in vizinhos[:3]]
    det = (v[0].x * (v[1].y * v[2].z - v[1].z * v[2].y)
           - v[0].y * (v[1].x * v[2].z - v[1].z * v[2].x)
           + v[0].z * (v[1].x * v[2].y - v[1].y * v[2].x))
    return 1 if det > 0 else -1


def referencia_3d(smiles):
    """Molecula 3D a partir do SMILES: L por construcao."""
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None, None
    mh = Chem.AddHs(m)
    if AllChem.EmbedMolecule(mh, AllChem.ETKDGv3()) != 0:
        return None, None
    AllChem.MMFFOptimizeMolecule(mh, mmffVariant="MMFF94s", maxIters=2000)
    m3 = Chem.RemoveHs(mh)
    Chem.AssignStereochemistryFrom3D(m3)
    return m3, m


def analisa(path, ref, forcado=None):
    base = os.path.splitext(os.path.basename(path))[0]
    alvo = forcado or nome_para_ref(base, ref)
    if alvo is None or alvo not in ref:
        return base, "?", "ERRO", "sem registro correspondente no arquivo de SMILES"

    m, erro = le_3d(path)
    if erro:
        return base, alvo, "ERRO", erro

    m3, _ = referencia_3d(ref[alvo]["smiles"])
    if m3 is None:
        return base, alvo, "ERRO", "nao foi possivel gerar a referencia 3D"

    # consulta com ligacoes genericas: o alvo nao tem ordens de ligacao
    p = Chem.AdjustQueryParameters()
    p.makeBondsGeneric = True
    p.makeDummiesQueries = True
    q = Chem.AdjustQueryProperties(Chem.Mol(m3), p)

    alvo_pesado = Chem.RemoveHs(Chem.Mol(m), sanitize=False)
    match = alvo_pesado.GetSubstructMatch(q)
    if not match:
        return base, alvo, "ERRO", "conectividade nao bate com a referencia"

    centros = Chem.FindMolChiralCenters(m3, useLegacyImplementation=False)
    conf_ref = m3.GetConformer()
    conf_alv = alvo_pesado.GetConformer()

    iguais, invertidos = [], []
    for pos, (idx, cip) in enumerate(centros, 1):
        viz = [a.GetIdx() for a in m3.GetAtomWithIdx(idx).GetNeighbors()]
        if len(viz) < 3:
            continue
        s_ref = sinal_volume(conf_ref, idx, viz)
        s_alv = sinal_volume(conf_alv, match[idx], [match[v] for v in viz])
        (iguais if s_ref == s_alv else invertidos).append((pos, cip))

    n = len(iguais) + len(invertidos)
    if invertidos:
        quais = ", ".join(f"centro {p} (era {c})" for p, c in invertidos)
        return base, alvo, "INVERTIDO", f"{len(invertidos)}/{n} invertido(s): {quais}"
    return base, alvo, "OK", f"{n}/{n} centros com a mesma configuracao da referencia (all-L)"


def main(args):
    forcado = None
    if "--ref" in args:
        i = args.index("--ref")
        forcado = args[i + 1]
        args = args[:i] + args[i + 2:]

    alvos = []
    for a in args:
        if os.path.isdir(a):
            alvos += [os.path.join(a, f) for f in sorted(os.listdir(a))
                      if os.path.splitext(f)[1].lower()
                      in (".pdb", ".ent", ".sdf", ".mol", ".mol2")]
        else:
            alvos.append(a)
    if not alvos:
        print(__doc__)
        return 1

    ref = carrega_referencia(REF)
    print(f"{'arquivo':16} {'ref':7} {'situacao':11} detalhe")
    print("-" * 88)
    n_ok = 0
    for p in alvos:
        base, alvo, sit, det = analisa(p, ref, forcado)
        n_ok += sit == "OK"
        print(f"{base:16} {alvo:7} {sit:11} {det}")
    print("-" * 88)
    print(f"{n_ok}/{len(alvos)} sem problema de estereoquimica")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
