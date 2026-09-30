"""Confere a estereoquimica de estruturas 3D contra os SMILES de referencia.

Serve para validar o que saiu de um gerador de 3D (por exemplo o tradutor web
que o ChimeraX usa no `open smiles:`) antes de usar as estruturas em docking:
se um centro quiral vier invertido, o residuo deixa de ser L e a molecula e
outra, com pose diferente no sitio. Minimizar depois nao corrige isso --
a minimizacao preserva a quiralidade que recebeu.

Uso:
    python3 scripts/checa_quiralidade_3D.py arquivo1.sdf arquivo2.pdb ...
    python3 scripts/checa_quiralidade_3D.py pasta_com_estruturas/

Formatos: .sdf/.mol (preferidos, tem ordem de ligacao explicita), .mol2, .pdb.
PDB e o menos confiavel: nao guarda ordem de ligacao, entao a percepcao pode
falhar em residuo fora do padrao. Se puder escolher, exporte SDF ou MOL2.

O nome do arquivo (sem extensao) e casado com o nome do registro em
Smiles_mcro_analogos.txt -- por exemplo mcr3.sdf casa com >mcr3.
"""
import os
import sys
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors

RDLogger.DisableLog('rdApp.*')

REF = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "..", "Smiles_mcro_analogos.txt")


def carrega_referencia(path):
    """Le o arquivo de SMILES: >nome / codigo / codigo* / SMILES."""
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
    saida = {}
    for nome, ls in ref.items():
        if len(ls) >= 3:
            saida[nome] = {"codigo": ls[1], "smiles": ls[2]}
    return saida


def le_3d(path):
    ext = os.path.splitext(path)[1].lower()
    if ext in (".sdf", ".mol"):
        m = next(iter(Chem.SDMolSupplier(path, removeHs=False)), None)
    elif ext == ".mol2":
        m = Chem.MolFromMol2File(path, removeHs=False)
    elif ext in (".pdb", ".ent"):
        m = Chem.MolFromPDBFile(path, removeHs=False)
    else:
        return None, f"extensao nao suportada ({ext})"
    if m is None:
        return None, "RDKit nao conseguiu ler o arquivo"
    if m.GetNumConformers() == 0:
        return None, "sem coordenadas 3D"
    return m, None


def centros(m):
    """CIP de cada centro quiral, deduzido das coordenadas 3D."""
    Chem.AssignStereochemistryFrom3D(m)
    return [(a.GetIdx(), a.GetPropsAsDict().get("_CIPCode", "?"))
            for a in m.GetAtoms() if a.HasProp("_CIPCode")]


def analisa(path, ref):
    nome = os.path.splitext(os.path.basename(path))[0]
    m, erro = le_3d(path)
    if erro:
        return nome, "ERRO", erro, None

    cips = centros(m)
    labs = [c for _, c in cips]
    n_s, n_r = labs.count("S"), labs.count("R")
    formula = rdMolDescriptors.CalcMolFormula(m)

    notas = []
    r = ref.get(nome)
    if r:
        mref = Chem.MolFromSmiles(r["smiles"])
        mref = Chem.AddHs(mref)
        f_ref = rdMolDescriptors.CalcMolFormula(mref)
        if formula.replace("+", "").replace("-", "") != f_ref:
            notas.append(f"formula difere da referencia ({formula} vs {f_ref})")
        # conectividade: compara o grafo sem estereo
        a = Chem.MolToSmiles(Chem.RemoveHs(Chem.MolFromSmiles(
            Chem.MolToSmiles(Chem.RemoveHs(m), isomericSmiles=False))))
        b = Chem.MolToSmiles(Chem.MolFromSmiles(
            Chem.MolToSmiles(Chem.MolFromSmiles(r["smiles"]),
                             isomericSmiles=False)))
        if a != b:
            notas.append("conectividade difere da referencia")
    else:
        notas.append("sem registro correspondente no arquivo de SMILES")

    if len(cips) != 6:
        notas.append(f"{len(cips)} centros quirais detectados (esperado 6)")
    if n_r:
        pos = [str(i + 1) for i, (_, c) in enumerate(cips) if c == "R"]
        notas.append(f"{n_r} centro(s) R = residuo D: centro(s) {', '.join(pos)}")

    ok = not notas
    status = "OK" if ok else ("ATENCAO" if n_r or len(cips) != 6 else "ver nota")
    return nome, status, "; ".join(notas) or "all-L, confere com a referencia", \
        f"{n_s}S/{n_r}R"


def main(args):
    alvos = []
    for a in args:
        if os.path.isdir(a):
            for f in sorted(os.listdir(a)):
                if os.path.splitext(f)[1].lower() in (
                        ".sdf", ".mol", ".mol2", ".pdb", ".ent"):
                    alvos.append(os.path.join(a, f))
        else:
            alvos.append(a)
    if not alvos:
        print(__doc__)
        return 1

    ref = carrega_referencia(REF)
    print(f"referencia: {len(ref)} registros\n")
    print(f"{'estrutura':14} {'status':9} {'centros':9} nota")
    print("-" * 78)
    n_prob = 0
    for p in alvos:
        nome, status, nota, cent = analisa(p, ref)
        n_prob += status != "OK"
        print(f"{nome:14} {status:9} {str(cent or '-'):9} {nota}")
    print("-" * 78)
    print(f"{len(alvos) - n_prob}/{len(alvos)} sem problema")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
