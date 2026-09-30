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

def nat(s):
    return Chem.MolFromSmiles(s, sanitize=False).GetNumAtoms()

# tons pastel: realce legivel sem competir com os atomos
COL = {2: (0.98, 0.78, 0.55),   # pos 2  F / V / L
       3: (0.66, 0.82, 0.95),   # pos 3  S / P
       4: (0.72, 0.89, 0.72),   # pos 4  P / L
       5: (0.87, 0.77, 0.95)}   # pos 5  E / K

mols, hA, hB, legends = [], [], [], []
for name, a, b, c, d in rows:
    blocks = [HEAD, B2[a], B3[b], B4[c], B5[d], TAIL]
    m = Chem.MolFromSmiles("".join(blocks))
    off, acol, res = 0, {}, {}
    for i, blk in enumerate(blocks, 1):
        n = nat(blk)
        if i in COL:
            for idx in range(off, off + n):
                acol[idx] = COL[i]
                res[idx] = i
        off += n
    bcol = {}
    for bd in m.GetBonds():
        i, j = bd.GetBeginAtomIdx(), bd.GetEndAtomIdx()
        if res.get(i) is not None and res.get(i) == res.get(j):
            bcol[bd.GetIdx()] = COL[res[i]]
    mols.append(m); hA.append(acol); hB.append(bcol)
    legends.append(f"{name}    Z-{a}-{b}-{c}-{d}-N*")

# Template = backbone completo (pGlu -> 5 residuos -> amida C-terminal), com as
# cadeias laterais reduzidas a Gly. Ancorar nele, e nao no MCS, evita que as
# moleculas saiam giradas/espelhadas umas em relacao as outras: o MCS de 18
# atomos e curto e casa em regioes diferentes de cada molecula.
TEMPLATE = "O=C1CCC(N1)C(=O)NCC(=O)NCC(=O)NCC(=O)NCC(=O)NCC(N)=O"

# CoordGen no lugar do algoritmo padrao: com o backbone fixo, o layout padrao
# encaixa as cadeias laterais sem espaco e sobrepoe atomos no plano. Medindo
# pares nao-ligados mais proximos que 0.70 x o comprimento medio de ligacao,
# o padrao da 60 colisoes nos 24 compostos e o CoordGen da 4.
rdDepictor.SetPreferCoordGen(True)

tpl = Chem.MolFromSmiles(TEMPLATE)
rdDepictor.Compute2DCoords(tpl)


def orienta_n_para_esquerda(m, match):
    """Poe o N-terminal (pGlu) a esquerda, convencao para peptideos.

    A checagem e feita na molecula ja desenhada, e nao no template: orientar o
    template nao se propaga de forma confiavel para o resultado final. `match`
    mapeia os atomos do template, entao match[0] e o O do pGlu e match[-1] e o
    O da amida C-terminal. Tem de ser rotacao de 180 graus e nao espelhamento,
    que inverteria a leitura da estereoquimica.
    """
    conf = m.GetConformer()
    if conf.GetAtomPosition(match[0]).x <= conf.GetAtomPosition(match[-1]).x:
        return False
    for i in range(m.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        p.x, p.y = -p.x, -p.y
        conf.SetAtomPosition(i, p)
    return True

n_ok = n_gir = 0
for m in mols:
    match = m.GetSubstructMatch(tpl)
    if match:
        rdDepictor.GenerateDepictionMatching2DStructure(m, tpl, acceptFailure=True)
        n_ok += 1
        n_gir += orienta_n_para_esquerda(m, match)
    else:
        rdDepictor.Compute2DCoords(m)
print(f"backbone-template ({tpl.GetNumAtoms()} atomos) casou em {n_ok}/{len(mols)}")
print(f"girados 180 graus para o N-terminal ficar a esquerda: {n_gir}/{len(mols)}")

def render(sub, path, ncols, w=640, h=520, svg=False):
    ms = [mols[i] for i in sub]
    A  = [hA[i] for i in sub]
    B  = [hB[i] for i in sub]
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

# blocos de 6 compostos, grade 3 x 2
CHUNK = 6
for k in range(0, len(rows), CHUNK):
    sub = list(range(k, min(k + CHUNK, len(rows))))
    nome = f"bloco{k//CHUNK + 1}_{rows[sub[0]][0]}-{rows[sub[-1]][0]}"
    render(sub, f"{OUT}/{nome}.png", 3)
    render(sub, f"{OUT}/{nome}.svg", 3, svg=True)

# um arquivo por composto, com realce contInuo
for i, (name, *_) in enumerate(rows):
    render([i], f"{OUT}/{name}.png", 1, 960, 760)
    render([i], f"{OUT}/{name}.svg", 1, 960, 760, svg=True)

print("pronto")
