# Dinâmica Molecular

## POSE DOCKING - L_mcro_TRPA1_dk2

### Parâmetros do sistema

- Campo de força: CHARMM36 (pasta `charmm36.ff` que veio no zip do CGenFF, usada para proteína **e** ligante)
- Caixa de água: 1.2 nm, modelo TIP3P
- Sais: NaCl [0,1 M] + neutralização
- Programa: GROMACS

---

## ETAPA 0 — Organizar a pasta de trabalho

Todos os comandos são rodados **de dentro desta mesma pasta**.

```
MD_dk2/
├── charmm36.ff/              # veio do zip do CGenFF (contém L-m_ffbonded.itp)
├── TRPA1_AF3_Y_dk2.pdb       # proteína (pose do docking, SEM o ligante)
├── L-mcro_dk2_gmx.pdb        # ligante (pose do docking)
├── L-mcro_dk2_gmx.top        # topologia do ligante gerada pelo CGenFF
└── ions.mdp  em.mdp  nvt.mdp  npt.mdp  md.mdp   # (conteúdo no final deste arquivo)
```

Conferir se o arquivo de parâmetros do ligante está dentro da pasta do campo de força:

```bash
ls charmm36.ff/L-m_ffbonded.itp
```

---

## ETAPA 1 — Protonação da proteína (pós-docking)

Programa PDB2PQR (https://server.poissonboltzmann.org/pdb2pqr), que gera os arquivos `.log` e `.pqr`.

Script Python para saber quais aminoácidos devem estar protonados (pH 7):

```bash
python3 Output-ProtonationPDB2PQR.py TRPA1_AF3_Y_dk2.log TRPA1_AF3_Y_dk2.pqr 7 > output_TRPA1_AF3_Y_dk2.txt
```

O `output_TRPA1_AF3_Y_dk2.txt` é usado para responder às perguntas do `pdb2gmx` na Etapa 2.

---

## ETAPA 2 — Topologia da proteína

```bash
gmx pdb2gmx -f TRPA1_AF3_Y_dk2.pdb -ff charmm36 -water tip3p -ignh -lys -his -asp -glu -ter -o protein.gro
```

- Se aparecer uma lista de campos de força, escolher o que diz **"from current directory"**.
- Responder LYS / HIS / ASP / GLU / terminais conforme o `output_TRPA1_AF3_Y_dk2.txt`.
- Arquivos gerados: `protein.gro`, `topol.top`, `posre.itp` (ou `topol_Protein_chain_X.itp` e `posre_Protein_chain_X.itp` se houver várias cadeias).

---

## ETAPA 3 — Topologia do ligante (CGenFF → GROMACS)

Topologia preparada no CGenFF (https://cgenff.com/) na saída GROMACS. Arquivos: `L-mcro_dk2_gmx.pdb`, `L-mcro_dk2_gmx.top` e `charmm36.ff/L-m_ffbonded.itp`.

> ⚠️ Antes de seguir: verificar as **penalidades** no arquivo `.str` do CGenFF (acima de ~10, revisar; acima de ~50, validar) e a **protonação do ligante** em pH 7 (a topologia atual tem carga total 0 e aminas neutras).

**3.1 — Extrair só a molécula do `.top` para um `.itp`** (do `[ moleculetype ]` até o fim de `[ dihedrals ]`):

```bash
sed -n '12,1186p' L-mcro_dk2_gmx.top > L-mcro_dk2_gmx.itp
```

**3.2 — Renomear a molécula e o resíduo para `LMC`** (sem hífen, igual no `.itp` e no `.pdb`):

```bash
sed -i 's/^Other /LMC   /; s/ L-m / LMC /g' L-mcro_dk2_gmx.itp
sed -i 's/ L-m / LMC /' L-mcro_dk2_gmx.pdb
```

**3.3 — Conferir:**

```bash
grep -A2 moleculetype L-mcro_dk2_gmx.itp     # deve mostrar: LMC   3
grep -c LMC L-mcro_dk2_gmx.pdb               # deve mostrar: 116
```

**3.4 — Converter o ligante para `.gro`:**

```bash
gmx editconf -f L-mcro_dk2_gmx.pdb -o lig.gro
```

**3.5 — Restrição de posição do ligante** (usada no NVT/NPT):

```bash
gmx make_ndx -f lig.gro -o index_lig.ndx
> 0 & ! a H*
> q
gmx genrestr -f lig.gro -n index_lig.ndx -o posre_lig.itp -fc 1000 1000 1000
```

No `genrestr`, escolher o grupo novo (`System_&_!H*`, átomos pesados).

---

## ETAPA 4 — Editar o `topol.top`

Abrir com `nano topol.top` e fazer **4 alterações**:

**4.1 — Logo depois do forcefield, incluir os parâmetros do ligante:**

```
; Include forcefield parameters
#include "./charmm36.ff/forcefield.itp"

; Include ligand parameters
#include "./charmm36.ff/L-m_ffbonded.itp"
```

**4.2 e 4.3 — Depois da proteína (depois do `posre.itp` / de todas as cadeias) e ANTES da água, incluir a topologia e a restrição do ligante:**

```
; Include Position restraint file
#ifdef POSRES
#include "posre.itp"
#endif

; Include ligand topology
#include "L-mcro_dk2_gmx.itp"

; Ligand position restraints
#ifdef POSRES_LIG
#include "posre_lig.itp"
#endif

; Include water topology
#include "./charmm36.ff/tip3p.itp"
```

**4.4 — No final, em `[ molecules ]`, adicionar o ligante logo após a proteína:**

```
[ molecules ]
; Compound        #mols
Protein_chain_A     1
LMC                 1
```

Salvar: `Ctrl + O`, `Enter`, `Ctrl + X`.

---

## ETAPA 5 — Montar o complexo proteína + ligante

```bash
head -n -1 protein.gro > complex.gro           # proteína sem a linha da caixa
sed -n '3,118p' lig.gro >> complex.gro          # 116 átomos do ligante
tail -n 1 protein.gro >> complex.gro            # linha da caixa
N=$(( $(sed -n 2p protein.gro) + 116 ))
sed -i "2s/.*/$N/" complex.gro                  # corrige o número total de átomos
```

Conferir visualmente (PyMOL/VMD) se o ligante está no sítio do docking:

```bash
gmx editconf -f complex.gro -o complex_check.pdb
```

---

## ETAPA 6 — Caixa (1.2 nm)

```bash
gmx editconf -f complex.gro -o box.gro -c -d 1.2 -bt dodecahedron
```

(Para caixa cúbica: `-bt cubic`.)

---

## ETAPA 7 — Solvatação (TIP3P)

```bash
gmx solvate -cp box.gro -cs spc216.gro -p topol.top -o solv.gro
```

---

## ETAPA 8 — Íons (NaCl 0,1 M + neutralização)

```bash
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 1
gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.1
```

Grupo a substituir: **SOL**.

---

## ETAPA 9 — Minimização de energia

```bash
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em
```

Conferir: `Fmax < 1000` e energia potencial negativa no final.

```bash
gmx energy -f em.edr -o potential.xvg      # escolher "Potential"
```

---

## ETAPA 10 — Grupos de índice (Protein_LMC)

```bash
gmx make_ndx -f em.gro -o index.ndx
> 1 | 13        # 1 = Protein; 13 = LMC (conferir o número do grupo LMC na lista)
> q
```

Isso cria o grupo `Protein_LMC`. `Water_and_ions` já existe por padrão.

---

## ETAPA 11 — Equilíbrio NVT (100 ps)

```bash
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
gmx mdrun -v -deffnm nvt
```

```bash
gmx energy -f nvt.edr -o temperature.xvg   # escolher "Temperature"
```

---

## ETAPA 12 — Equilíbrio NPT (100 ps)

```bash
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -n index.ndx -o npt.tpr
gmx mdrun -v -deffnm npt
```

```bash
gmx energy -f npt.edr -o pressure.xvg      # escolher "Pressure"
gmx energy -f npt.edr -o density.xvg       # escolher "Density"
```

---

## ETAPA 13 — Produção (100 ns)

```bash
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md.tpr
gmx mdrun -v -deffnm md
```

Se a simulação parar, continuar de onde parou:

```bash
gmx mdrun -v -deffnm md -cpi md.cpt
```

---

## ETAPA 14 — Pós-processamento básico

Centralizar e corrigir a periodicidade:

```bash
gmx trjconv -s md.tpr -f md.xtc -n index.ndx -o md_center.xtc -center -pbc mol -ur compact
# centralizar: Protein_LMC | saída: System
```

RMSD da proteína e do ligante:

```bash
gmx rms -s md.tpr -f md_center.xtc -n index.ndx -o rmsd_protein.xvg -tu ns   # Backbone / Backbone
gmx rms -s md.tpr -f md_center.xtc -n index.ndx -o rmsd_lig.xvg -tu ns       # Backbone / LMC
```

---

## ARQUIVOS .mdp

### ions.mdp

```
integrator      = steep
emtol           = 1000.0
emstep          = 0.01
nsteps          = 50000
nstlist         = 1
cutoff-scheme   = Verlet
coulombtype     = cutoff
rcoulomb        = 1.2
rvdw            = 1.2
pbc             = xyz
```

### em.mdp

```
integrator      = steep
emtol           = 1000.0
emstep          = 0.01
nsteps          = 50000
nstlist         = 10
cutoff-scheme   = Verlet
coulombtype     = PME
rcoulomb        = 1.2
vdwtype         = Cut-off
vdw-modifier    = Force-switch
rvdw-switch     = 1.0
rvdw            = 1.2
DispCorr        = no
pbc             = xyz
```

### nvt.mdp

```
define                  = -DPOSRES -DPOSRES_LIG
integrator              = md
nsteps                  = 50000         ; 2 fs * 50000 = 100 ps
dt                      = 0.002
nstxout-compressed      = 5000
nstenergy               = 500
nstlog                  = 500
continuation            = no
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 1
lincs_order             = 4
cutoff-scheme           = Verlet
nstlist                 = 20
coulombtype             = PME
pme_order               = 4
fourierspacing          = 0.16
rcoulomb                = 1.2
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2
DispCorr                = no
tcoupl                  = V-rescale
tc-grps                 = Protein_LMC Water_and_ions
tau_t                   = 0.1   0.1
ref_t                   = 300   300
pcoupl                  = no
pbc                     = xyz
gen_vel                 = yes
gen_temp                = 300
gen_seed                = -1
```

### npt.mdp

```
define                  = -DPOSRES -DPOSRES_LIG
integrator              = md
nsteps                  = 50000         ; 100 ps
dt                      = 0.002
nstxout-compressed      = 5000
nstenergy               = 500
nstlog                  = 500
continuation            = yes
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 1
lincs_order             = 4
cutoff-scheme           = Verlet
nstlist                 = 20
coulombtype             = PME
pme_order               = 4
fourierspacing          = 0.16
rcoulomb                = 1.2
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2
DispCorr                = no
tcoupl                  = V-rescale
tc-grps                 = Protein_LMC Water_and_ions
tau_t                   = 0.1   0.1
ref_t                   = 300   300
pcoupl                  = C-rescale
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = 1.0
compressibility         = 4.5e-5
refcoord_scaling        = com
pbc                     = xyz
gen_vel                 = no
```

### md.mdp

```
integrator              = md
nsteps                  = 50000000      ; 2 fs * 50.000.000 = 100 ns
dt                      = 0.002
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstxout-compressed      = 5000          ; salva a cada 10 ps
compressed-x-grps       = System
nstenergy               = 5000
nstlog                  = 5000
continuation            = yes
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 1
lincs_order             = 4
cutoff-scheme           = Verlet
nstlist                 = 20
coulombtype             = PME
pme_order               = 4
fourierspacing          = 0.16
rcoulomb                = 1.2
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2
DispCorr                = no
tcoupl                  = V-rescale
tc-grps                 = Protein_LMC Water_and_ions
tau_t                   = 0.1   0.1
ref_t                   = 300   300
pcoupl                  = C-rescale
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = 1.0
compressibility         = 4.5e-5
pbc                     = xyz
gen_vel                 = no
```

> Observações:
> - `C-rescale` exige GROMACS 2021 ou mais novo. Em versões antigas, use `Berendsen` no NPT e `Parrinello-Rahman` na produção.
> - Temperatura em 300 K; para condição fisiológica, trocar `ref_t` e `gen_temp` para 310.
> - TRPA1 é canal de membrana: se o modelo incluir a região transmembrana, considerar montar o sistema em bicamada lipídica (CHARMM-GUI Membrane Builder).
