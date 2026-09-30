"""Gera as imagens 2D dos 24 analogos, em blocos de 6 e individuais.

Layout: alinhamento por MCS, com o algoritmo de depicao padrao do RDKit.
Medindo pares de atomos nao-ligados mais proximos que 0.70 x o comprimento
medio de ligacao, essa combinacao da 0 colisoes nos 24 compostos. As
alternativas testadas sao piores: CoordGen livre da 8, template do backbone
com CoordGen da 4, e template do backbone sem CoordGen da 60 -- prender o
backbone tira o espaco de que as cadeias laterais precisam.

O custo do MCS e que as moleculas nao ficam todas na mesma orientacao entre
os paineis, ao contrario do que o template do backbone garantia.

Use scripts/mede_colisoes_2D.py para reavaliar se mudar o layout.
"""
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFMCS, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

RDLogger.DisableLog('rdApp.*')

OUT = "/home/user/Cloude/imagens_2D"

HEAD = "C1CC(=O)N[C@@H]1C(=O)"
TAIL = "N[C@@H](CC(N)=O)C(N)=O"
B2 = {"F": "N[C@@H](CC2=CC=CC=C2)C(=O)",
      "V": "N[C@@H](C(C)C)C(=O)",
      "L": "N[C@@H](CC(C)C)C(=O)"}
B3 = {"S": "N[C@@H](CO)C(=O)", "P": "N1CCC[C@H]1C(=O)"}
B4 = {"P": "N3CCC[C@H]3C(=O)", "L": "N[C@@H](CC(C)C)C(=O)"}
B5 = {"E": "N[C@@H](CCC(=O)O)C(=O)", "K": "N[C@@H](CCCCN)C(=O)"}

rows = [("mcro","F","S","P","E"),("mcr1","F","P","P","E"),("mcr2","F","S","P","K"),
        ("mcr3","F","S","L","E"),("mcr4","F","P","P","K"),("mcr5","F","P","L","E"),
        ("mcr6","F","S","L","K"),("mcr7","F","P","L","K"),
        ("mcr8","V","S","P","E"),("mcr9","L","S","P","E"),("mcr10","V","P","P","E"),
        ("mcr11","L","P","P","E"),("mcr12","V","S","P","K"),("mcr13","L","S","P","K"),
        ("mcr14","V","S","L","E"),("mcr15","L","S","L","E"),("mcr16","V","P","P","K"),
        ("mcr17","L","P","P","K"),("mcr18","V","P","L","E"),("mcr19","L","P","L","E"),
        ("mcr20","V","S","L","K"),("mcr21","L","S","L","K"),("mcr22","V","P","L","K"),
        ("mcr23","L","P","L","K")]

# tons pastel: realce legivel sem competir com os rotulos dos atomos
COL = {2: (0.98, 0.78, 0.55),   # pos 2  Phe / Val / Leu
       3: (0.66, 0.82, 0.95),   # pos 3  Ser / Pro
       4: (0.72, 0.89, 0.72),   # pos 4  Pro / Leu
       5: (0.87, 0.77, 0.95)}   # pos 5  Glu / Lys


def nat(smi):
    return Chem.MolFromSmiles(smi, sanitize=False).GetNumAtoms()


mols, hA, hB, legends = [], [], [], []
for name, a, b, c, d in rows:
    blocks = [HEAD, B2[a], B3[b], B4[c], B5[d], TAIL]
    m = Chem.MolFromSmiles("".join(blocks))
    # os indices seguem a ordem do SMILES, entao as faixas por residuo saem
    # do numero de atomos de cada bloco
    off, res = 0, {}
    for i, blk in enumerate(blocks, 1):
        n = nat(blk)
        if i in COL:
            for idx in range(off, off + n):
                res[idx] = i
        off += n
    acol = {i: COL[r] for i, r in res.items()}
    bcol = {}
    for bd in m.GetBonds():
        i, j = bd.GetBeginAtomIdx(), bd.GetEndAtomIdx()
        if res.get(i) is not None and res.get(i) == res.get(j):
            bcol[bd.GetIdx()] = COL[res[i]]
    mols.append(m); hA.append(acol); hB.append(bcol)
    legends.append(f"{name}    Z-{a}-{b}-{c}-{d}-N*")

mcs = rdFMCS.FindMCS(mols, ringMatchesRingOnly=True, timeout=120)
patt = Chem.MolFromSmarts(mcs.smartsString)
ref = Chem.Mol(mols[0])
rdDepictor.Compute2DCoords(ref)
for m in mols:
    try:
        rdDepictor.GenerateDepictionMatching2DStructure(m, ref, refPatt=patt)
    except Exception:
        rdDepictor.Compute2DCoords(m)
print(f"MCS comum: {mcs.numAtoms} atomos / {mcs.numBonds} ligacoes")


def render(sub, path, ncols, w=640, h=520, svg=False):
    ms = [mols[i] for i in sub]
    A = [hA[i] for i in sub]
    B = [hB[i] for i in sub]
    lg = [legends[i] for i in sub]
    nrows = -(-len(ms) // ncols)
    drv = rdMolDraw2D.MolDraw2DSVG if svg else rdMolDraw2D.MolDraw2DCairo
    d = drv(w * ncols, h * nrows, w, h)
    o = d.drawOptions()
    o.addStereoAnnotation = False
    o.bondLineWidth = 2
    o.legendFontSize = 30
    o.fillHighlights = True
    o.highlightBondWidthMultiplier = 20
    o.clearBackground = True
    d.DrawMolecules(ms, legends=lg,
                    highlightAtoms=[list(x) for x in A], highlightAtomColors=A,
                    highlightBonds=[list(x) for x in B], highlightBondColors=B)
    d.FinishDrawing()
    t = d.GetDrawingText()
    open(path, "w" if svg else "wb").write(t)
    print("  ", path)


CHUNK = 6
for k in range(0, len(rows), CHUNK):
    sub = list(range(k, min(k + CHUNK, len(rows))))
    nome = f"bloco{k // CHUNK + 1}_{rows[sub[0]][0]}-{rows[sub[-1]][0]}"
    render(sub, f"{OUT}/{nome}.png", 3)
    render(sub, f"{OUT}/{nome}.svg", 3, svg=True)

for i, (name, *_) in enumerate(rows):
    render([i], f"{OUT}/{name}.png", 1, 960, 760)
    render([i], f"{OUT}/{name}.svg", 1, 960, 760, svg=True)

print("pronto")
