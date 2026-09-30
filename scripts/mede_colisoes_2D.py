"""Compara layouts 2D medindo colisoes de atomos no plano."""
import math
from rdkit import Chem, RDLogger
from rdkit.Chem import rdDepictor
RDLogger.DisableLog('rdApp.*')

HEAD = "C1CC(=O)N[C@@H]1C(=O)"
TAIL = "N[C@@H](CC(N)=O)C(N)=O"
B2 = {"F": "N[C@@H](CC2=CC=CC=C2)C(=O)", "V": "N[C@@H](C(C)C)C(=O)",
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

TEMPLATE = "O=C1CCC(N1)C(=O)NCC(=O)NCC(=O)NCC(=O)NCC(=O)NCC(N)=O"


def build(name, a, b, c, d):
    return Chem.MolFromSmiles(HEAD + B2[a] + B3[b] + B4[c] + B5[d] + TAIL)


def collisions(m, frac=0.70):
    """Pares de atomos nao-ligados mais proximos que frac * comprimento medio
    de ligacao. Mede o quanto a estrutura esta sobreposta no plano."""
    conf = m.GetConformer()
    pts = [conf.GetAtomPosition(i) for i in range(m.GetNumAtoms())]
    bl = []
    for b in m.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        bl.append(math.hypot(pts[i].x - pts[j].x, pts[i].y - pts[j].y))
    mean_bl = sum(bl) / len(bl)
    cut = frac * mean_bl
    bonded = {(min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
               max(b.GetBeginAtomIdx(), b.GetEndAtomIdx())) for b in m.GetBonds()}
    n = 0
    for i in range(m.GetNumAtoms()):
        for j in range(i + 1, m.GetNumAtoms()):
            if (i, j) in bonded:
                continue
            if math.hypot(pts[i].x - pts[j].x, pts[i].y - pts[j].y) < cut:
                n += 1
    return n


def make_template(coordgen):
    rdDepictor.SetPreferCoordGen(coordgen)
    t = Chem.MolFromSmiles(TEMPLATE)
    rdDepictor.Compute2DCoords(t)
    conf = t.GetConformer()
    for i in range(t.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        p.x, p.y = -p.x, -p.y
        conf.SetAtomPosition(i, p)
    return t


def run(coordgen, constrained):
    rdDepictor.SetPreferCoordGen(coordgen)
    tpl = make_template(coordgen) if constrained else None
    rdDepictor.SetPreferCoordGen(coordgen)
    out = {}
    for name, a, b, c, d in rows:
        m = build(name, a, b, c, d)
        if constrained and m.HasSubstructMatch(tpl):
            rdDepictor.GenerateDepictionMatching2DStructure(m, tpl, acceptFailure=True)
        else:
            rdDepictor.Compute2DCoords(m)
        out[name] = collisions(m)
    return out


def run_mcs(coordgen):
    """Alinhamento por MCS -- o layout usado em gera_imagens_2D.py."""
    from rdkit.Chem import rdFMCS
    rdDepictor.SetPreferCoordGen(coordgen)
    ms = [build(*r) for r in rows]
    mcs = rdFMCS.FindMCS(ms, ringMatchesRingOnly=True, timeout=120)
    patt = Chem.MolFromSmarts(mcs.smartsString)
    ref = Chem.Mol(ms[0])
    rdDepictor.Compute2DCoords(ref)
    for m in ms:
        try:
            rdDepictor.GenerateDepictionMatching2DStructure(m, ref, refPatt=patt)
        except Exception:
            rdDepictor.Compute2DCoords(m)
    return {r[0]: collisions(m) for r, m in zip(rows, ms)}


variants = {
    "MCS, sem coordgen  (EM USO)    ": run_mcs(False),
    "MCS + coordgen                 ": run_mcs(True),
    "template backbone, sem coordgen": run(False, True),
    "template backbone + coordgen   ": run(True, True),
    "coordgen livre (sem template)  ": run(True, False),
}

print(f"{'variante':34} {'total':>6} {'piores'}")
for k, v in variants.items():
    tot = sum(v.values())
    piores = sorted(v.items(), key=lambda x: -x[1])[:5]
    piores = ", ".join(f"{n}:{c}" for n, c in piores if c > 0) or "nenhum"
    print(f"{k} {tot:6} {piores}")

print("\ndetalhe por composto (colisoes):")
print(f"{'composto':10} {'MCS':>6} {'tpl':>7} {'cg livre':>9}")
for name, *_ in rows:
    a = variants["MCS, sem coordgen  (EM USO)    "][name]
    b = variants["template backbone, sem coordgen"][name]
    c = variants["coordgen livre (sem template)  "][name]
    flag = "  <--" if a > 0 and b == 0 else ""
    print(f"{name:10} {a:6} {b:7} {c:9}{flag}")
