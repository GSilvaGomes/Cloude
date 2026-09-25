# Dinâmica Molecular

## POSE DOCKING - L_mcro_TRPA1_dk2

### Parâmetros do sistema

- Campo de força: CHARMM36 (pasta `charmm36.ff` que veio no zip do CGenFF, usada para proteína **e** ligante)
- Caixa de água: 1.2 nm, modelo TIP3P
- Sais: NaCl [0,1 M] + neutralização
- Programa: GROMACS

---

## ETAPA 0 — Organizar a pasta de trabalho

Todos os comandos são rodados **de dentro da pasta do sistema** (`MD_dk2/`). Os arquivos `.mdp` ficam na pasta de cima (`../`).

```
pasta_mdp/
├── ions.mdp  minim.mdp  nvt.mdp  npt.mdp  md.mdp   # (conteúdo no final deste arquivo)
└── MD_dk2/
    ├── charmm36.ff/              # veio do zip do CGenFF (contém L-m_ffbonded.itp)
    ├── TRPA1_AF3_Y_dk2.pdb       # proteína (pose do docking, SEM o ligante)
    ├── L-mcro_dk2_gmx.pdb        # ligante (pose do docking)
    └── L-mcro_dk2_gmx.top        # topologia do ligante gerada pelo CGenFF
```

Conferir se o arquivo de parâmetros do ligante está dentro da pasta do campo de força:

```bash
ls charmm36.ff/L-m_ffbonded.itp
```

> No Vital o executável é `gmx` (GROMACS 2022, com GPU). Usar sempre o mesmo em todas as etapas.

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
gmx pdb2gmx -f TRPA1_AF3_Y_dk2.pdb -ff charmm36 -water tip3p -ignh -lys -his -asp -glu -ter -o 0-protein.gro
```

- Se aparecer uma lista de campos de força, escolher o que diz **"from current directory"**.
- Responder LYS / HIS / ASP / GLU / terminais conforme o `output_TRPA1_AF3_Y_dk2.txt`.
- Arquivos gerados: `0-protein.gro`, `topol.top`, `posre.itp` (ou `topol_Protein_chain_X.itp` e `posre_Protein_chain_X.itp` se houver várias cadeias).

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

**3.5 — Restrição de posição do ligante** (usada no NVT/NPT; fazer logo após a 3.4 e antes da Etapa 4):

```bash
gmx make_ndx -f lig.gro -o index_lig.ndx
> 0 & ! a H*        # cria o grupo 3 "System_&_!H*" (só átomos pesados)
> q
gmx genrestr -f lig.gro -n index_lig.ndx -o posre_lig.itp -fc 1000 1000 1000
# escolher o grupo 3 (System_&_!H*)
```

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

**4.4 — No final, em `[ molecules ]`, adicionar o ligante logo após a proteína (depois da última cadeia):**

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
head -n -1 0-protein.gro > 0-complex.gro        # proteína sem a linha da caixa
sed -n '3,118p' lig.gro >> 0-complex.gro         # 116 átomos do ligante
tail -n 1 0-protein.gro >> 0-complex.gro         # linha da caixa
N=$(( $(sed -n 2p 0-protein.gro) + 116 ))
sed -i "2s/.*/$N/" 0-complex.gro                 # corrige o número total de átomos
```

Conferência visual (opcional) — abrir no PyMOL/VMD e ver se o ligante está no sítio do docking, inteiro e sem sobreposição com a proteína:

```bash
gmx editconf -f 0-complex.gro -o 0-complex_check.pdb
```

---

## ETAPA 6 — Caixa (1.2 nm)

⚠️ Usar o **complexo** (`0-complex.gro`), não o `0-protein.gro` — senão o ligante fica fora do sistema.

```bash
gmx editconf -f 0-complex.gro -o 1-box.gro -bt triclinic -d 1.2 -c
```

(Alternativa com menos água e mais rápida: `-bt dodecahedron`.)

---

## ETAPA 7 — Solvatação (TIP3P)

```bash
gmx solvate -cp 1-box.gro -cs spc216.gro -p topol.top -o 2-solvate.gro
```

---

## ETAPA 8 — Íons (NaCl 0,1 M + neutralização)

```bash
gmx grompp -f ../ions.mdp -c 2-solvate.gro -p topol.top -o 3-ions.tpr -po 3-ions-out.mdp
gmx genion -s 3-ions.tpr -o 4-solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.1
# escolher o grupo SOL
```

Só se o `grompp` parar por aviso de carga total diferente de zero, acrescentar `-maxwarn 1`.

---

## ETAPA 9 — Grupos de índice (Protein_LMC)

Criar o `index.ndx` **antes da minimização**, porque todos os `grompp` seguintes usam `-n index.ndx`:

```bash
gmx make_ndx -f 4-solv_ions.gro -o index.ndx
> 1 | 13        # 1 = Protein; 13 = LMC (conferir o número do grupo LMC na lista)
> q
```

Isso cria o grupo `Protein_LMC`. `Water_and_ions` já existe por padrão.

---

## ETAPA 10 — Minimização de energia

```bash
gmx grompp -f ../minim.mdp -c 4-solv_ions.gro -p topol.top -n index.ndx -o 5-minim.tpr -po 5-minim-out.mdp
gmx mdrun -v -deffnm 5-minim -ntmpi 1 > 5-minim_run.log 2>&1
```

Conferir: `Fmax < 1000` e energia potencial negativa no final.

```bash
gmx energy -f 5-minim.edr -o 5-potential.xvg      # escolher "Potential"
```

---

## ETAPA 11 — Equilíbrio NVT (100 ps)

```bash
gmx grompp -f ../nvt.mdp -c 5-minim.gro -r 5-minim.gro -p topol.top -n index.ndx -o 6-nvt.tpr -po 6-nvt-out.mdp
gmx mdrun -v -deffnm 6-nvt -ntmpi 1 > 6-nvt_run.log 2>&1
```

```bash
gmx energy -f 6-nvt.edr -o 6-temperature.xvg      # escolher "Temperature"
```

---

## ETAPA 12 — Equilíbrio NPT (100 ps)

```bash
gmx grompp -f ../npt.mdp -c 6-nvt.gro -r 6-nvt.gro -t 6-nvt.cpt -p topol.top -n index.ndx -o 7-npt.tpr -po 7-npt-out.mdp
gmx mdrun -v -deffnm 7-npt -ntmpi 1 > 7-npt_run.log 2>&1
```

```bash
gmx energy -f 7-npt.edr -o 7-pressure.xvg         # escolher "Pressure"
gmx energy -f 7-npt.edr -o 7-density.xvg          # escolher "Density"
```

---

## ETAPA 13 — Produção (100 ns)

```bash
gmx grompp -f ../md.mdp -c 7-npt.gro -t 7-npt.cpt -p topol.top -n index.ndx -o 8-md.tpr -po 8-md-out.mdp
gmx mdrun -v -deffnm 8-md -ntmpi 1 > 8-md_run.log 2>&1
```

Se a simulação parar, continuar de onde parou:

```bash
gmx mdrun -v -deffnm 8-md -cpi 8-md.cpt -ntmpi 1 >> 8-md_run.log 2>&1
```

---

## ETAPA 14 — Pós-processamento básico

Centralizar e corrigir a periodicidade:

```bash
gmx trjconv -s 8-md.tpr -f 8-md.xtc -n index.ndx -o 9-md_center.xtc -center -pbc mol -ur compact
# centralizar: Protein_LMC | saída: System
```

RMSD da proteína e do ligante:

```bash
gmx rms -s 8-md.tpr -f 9-md_center.xtc -n index.ndx -o 9-rmsd_protein.xvg -tu ns   # Backbone / Backbone
gmx rms -s 8-md.tpr -f 9-md_center.xtc -n index.ndx -o 9-rmsd_lig.xvg -tu ns       # Backbone / LMC
```

---

## ETAPA 15 — Rodar da minimização à produção no Slurm (servidor Vital)

As Etapas 0 a 9 são feitas **à mão no terminal** porque têm perguntas interativas (`pdb2gmx`, `genion`, `make_ndx`). Da minimização (Etapa 10) até a produção (Etapa 13) não há perguntas, então tudo roda de uma vez em um **job do Slurm**.

Antes de submeter, conferir se estes arquivos estão na pasta do sistema: `4-solv_ions.gro`, `topol.top`, `index.ndx`, `posre*.itp`, `L-mcro_dk2_gmx.itp`, `charmm36.ff/` — e os `.mdp` em `../`.

### Passo a passo no terminal do Ubuntu

**1. Abrir o terminal:** `Ctrl + Alt + T`.

**2. Entrar no servidor Vital por SSH** (trocar pelo seu usuário e o endereço do Vital):

```bash
ssh seu_usuario@endereco_do_vital
```

Digitar a senha (ela não aparece enquanto você digita; é normal) e apertar `Enter`.

**3. Ir até a pasta do sistema:**

```bash
cd ~/caminho/pasta_mdp/MD_dk2
pwd          # confirma em que pasta você está
ls           # confirma se 4-solv_ions.gro, topol.top, index.ndx... estão aqui
ls ../*.mdp  # confirma se os .mdp estão na pasta de cima
```

**4. Criar o arquivo do job com o nano:**

```bash
nano md_dk2.job
```

Abre uma tela de edição vazia.

**5. Colar o script** (conteúdo do item 15.2 abaixo): copiar o texto e, no terminal, colar com `Ctrl + Shift + V` (ou botão direito → Colar). No terminal, `Ctrl + V` sozinho não cola.

**6. Salvar e sair do nano:**
- `Ctrl + O` → aparece `File Name to Write: md_dk2.job` → apertar `Enter`
- `Ctrl + X` para sair

**7. Conferir se o arquivo ficou certo:**

```bash
cat md_dk2.job
```

Se o script foi copiado de um arquivo do Windows, corrigir as quebras de linha (senão o Slurm dá erro estranho):

```bash
sed -i 's/\r$//' md_dk2.job
```

**8. Submeter e acompanhar** (itens 15.3 e 15.4 abaixo):

```bash
sbatch md_dk2.job
squeue -u $USER
```

**9. Pode fechar o terminal.** O job continua rodando no servidor. Para voltar e ver o andamento, repetir os passos 1–3 e usar `cat md_dk2.out` / `tail -f 8-md.log`.

> **Alternativa:** criar o `md_dk2.job` no seu computador (com `nano md_dk2.job` no terminal local) e enviar para o Vital:
> ```bash
> scp md_dk2.job seu_usuario@endereco_do_vital:~/caminho/pasta_mdp/MD_dk2/
> ```

**15.1 — Criar o arquivo do job** (dentro da pasta do sistema):

```bash
nano md_dk2.job
```

**15.2 — Colar o conteúdo abaixo:**

```bash
#!/bin/bash -l

#SBATCH --job-name=md_dk2
#SBATCH --output=md_dk2.out
#SBATCH --error=md_dk2.err
#SBATCH --mem=64G
#SBATCH --time=168:00:00
#SBATCH --cpus-per-task=16

# Para o job se algum comando der erro (não segue para a próxima etapa)
set -e

# Entra na pasta de onde o job foi submetido
cd $SLURM_SUBMIT_DIR

# Se no Vital o GROMACS for carregado por módulo, descomentar e ajustar:
# module load gromacs

NT=${SLURM_CPUS_PER_TASK:-16}

# ---------- Verificações antes de começar ----------
# O ligante (116 átomos) precisa estar no grupo de temperatura junto com a proteína
grep -q "\[ Protein_LMC \]" index.ndx || { echo "ERRO: grupo Protein_LMC não existe no index.ndx (refazer Etapa 9)"; exit 1; }
for f in ../nvt.mdp ../npt.mdp ../md.mdp; do
  grep -Eq "^tc[-_]grps.*Protein_LMC" $f || { echo "ERRO: tc-grps em $f não usa Protein_LMC"; exit 1; }
  if grep -Eq "^comm[-_]grps" $f && ! grep -Eq "^comm[-_]grps.*(System|Protein_LMC)" $f; then
    echo "ERRO: comm-grps em $f não inclui o ligante (usar comm-grps = System)"; exit 1
  fi
done

echo "Início: $(date)"

# ---------- Minimização ----------
gmx grompp -f ../minim.mdp -c 4-solv_ions.gro -p topol.top -n index.ndx -o 5-minim.tpr -po 5-minim-out.mdp
gmx mdrun -v -deffnm 5-minim -ntmpi 1 -ntomp $NT > 5-minim_run.log 2>&1
echo "Minimização OK: $(date)"

# ---------- Equilíbrio NVT ----------
gmx grompp -f ../nvt.mdp -c 5-minim.gro -r 5-minim.gro -p topol.top -n index.ndx -o 6-nvt.tpr -po 6-nvt-out.mdp
gmx mdrun -v -deffnm 6-nvt -ntmpi 1 -ntomp $NT > 6-nvt_run.log 2>&1
echo "NVT OK: $(date)"

# ---------- Equilíbrio NPT ----------
gmx grompp -f ../npt.mdp -c 6-nvt.gro -r 6-nvt.gro -t 6-nvt.cpt -p topol.top -n index.ndx -o 7-npt.tpr -po 7-npt-out.mdp
gmx mdrun -v -deffnm 7-npt -ntmpi 1 -ntomp $NT > 7-npt_run.log 2>&1
echo "NPT OK: $(date)"

# ---------- Produção ----------
gmx grompp -f ../md.mdp -c 7-npt.gro -t 7-npt.cpt -p topol.top -n index.ndx -o 8-md.tpr -po 8-md-out.mdp
gmx mdrun -v -deffnm 8-md -ntmpi 1 -ntomp $NT > 8-md_run.log 2>&1
echo "Produção OK: $(date)"
```

Salvar e sair: `Ctrl + O`, `Enter`, `Ctrl + X`.

> - `--mem`: o tutorial do Vital recomenda bastante memória (512G), mas o GROMACS usa pouca; 64G é mais que suficiente e o job tende a entrar na fila mais rápido.
> - `--time` e `--cpus-per-task`: ajustar conforme os limites do Vital. O tempo precisa cobrir os 100 ns (ver quantos ns/dia aparecem no final do `7-npt.log` para estimar).
> - O `-v` no `mdrun` escreve a contagem de passos e a **previsão de término** no `_run.log` de cada etapa (ver 15.4).

**15.3 — Submeter o job:**

```bash
sbatch md_dk2.job
```

Aparece `Submitted batch job <número>`.

**15.4 — Acompanhar:**

```bash
squeue -u $USER          # ST: R = rodando | CD = terminou | F = falhou
cat md_dk2.out           # mostra quais etapas já terminaram (Minimização OK, NVT OK...)
cat md_dk2.err           # erros, se algo falhar
tail -f 8-md.log         # progresso da produção em tempo real (Ctrl + C sai)
```

Ver a contagem e a previsão de término da etapa que está rodando:

```bash
for f in 5-minim 6-nvt 7-npt 8-md; do [ -f ${f}_run.log ] && echo "$f: $(tr '\r' '\n' < ${f}_run.log | grep -E 'will finish|Step=|^step' | tail -1)"; done
```

Aparece algo como `8-md: step 2500000, will finish Sat Sep 27 14:32:10 2026`.
scancel <número_do_job>  # cancelar o job, se precisar
```

**15.5 — Se a produção parar no meio** (ex.: estourou o `--time`), criar um job de continuação `md_dk2_cont.job` com o mesmo cabeçalho `#SBATCH` (trocando `--job-name`, `--output` e `--error` para `md_dk2_cont`) e só estas linhas de execução:

```bash
cd $SLURM_SUBMIT_DIR
gmx mdrun -v -deffnm 8-md -cpi 8-md.cpt -ntmpi 1 -ntomp $SLURM_CPUS_PER_TASK >> 8-md_run.log 2>&1
```

E submeter com `sbatch md_dk2_cont.job`.

**15.6 — Rodar direto no seu Ubuntu (sem Slurm)**

Se for rodar num computador próprio com Ubuntu (sem `sbatch`), usar o mesmo script como `.sh`:

```bash
cp md_dk2.job md_dk2.sh
chmod +x md_dk2.sh                        # dá permissão para executar
nohup ./md_dk2.sh > md_dk2.out 2> md_dk2.err &
```

- As linhas `#SBATCH` são ignoradas (viram comentário).
- Trocar `NT=${SLURM_CPUS_PER_TASK:-16}` pelo número de núcleos do computador (ver com `nproc`), e `cd $SLURM_SUBMIT_DIR` por `cd "$(dirname "$0")"`.
- Se no seu computador o executável for `gmx_mpi` (e não `gmx`), trocar em todo o script e tirar o `-ntmpi 1` (que só funciona com `gmx`): `sed -i 's/gmx /gmx_mpi /g; s/ -ntmpi 1//' md_dk2.sh`.
- O `nohup ... &` deixa rodando mesmo se fechar o terminal. Acompanhar com `tail -f 8-md.log`; ver se ainda está rodando com `ps aux | grep mdrun`; parar com `kill <PID>`.

---

## PROBLEMAS ENCONTRADOS E SOLUÇÕES

### 1. Job falhou (`F`, `NonZeroExitCode`) em poucos segundos

Ver onde parou (na pasta do sistema):

```bash
cat md_dk2.out              # até qual etapa chegou
tail -40 md_dk2.err         # erros do grompp / command not found
tail -40 5-minim_run.log    # erros do mdrun (trocar pelo log da etapa que falhou)
```

### 2. `gmx` ou `gmx_mpi`?

No Vital o executável é **`gmx`** (`/usr/local/src/gromacs/gromacs/bin/gmx`, versão 2022). Usar `gmx` em todo o script:

```bash
sed -i 's/gmx_mpi/gmx/g' md_dk2.job
```

### 3. Erro no `mdrun`: GPU e `-ntomp` sem `-ntmpi`

```
Fatal error:
When using GPUs, setting the number of OpenMP threads without specifying the
number of ranks can lead to conflicting demands. Please specify the number of
thread-MPI ranks as well (option -ntmpi).
```

**Causa:** o nó do Vital tem GPU; com `-ntomp` o GROMACS exige também `-ntmpi`.

**Solução:** acrescentar `-ntmpi 1` em todos os `mdrun` (1 processo usando a GPU + as threads de CPU):

```bash
sed -i 's/ -ntomp / -ntmpi 1 -ntomp /' md_dk2.job
grep -n "mdrun" md_dk2.job
```

> `-ntmpi` só funciona com `gmx` (não com `gmx_mpi`).

### 4. Erro no `grompp` do NVT: 116 átomos fora dos grupos de temperatura

```
Fatal error:
116 atoms are not part of any of the T-Coupling groups
```

**Causa:** os 116 átomos são o ligante (LMC). O `tc-grps` dos `.mdp` não incluía o ligante (ex.: `Protein Water_and_ions`).

**Solução:**

a) Conferir os grupos usados e os que existem:

```bash
grep -n "tc-grps\|tc_grps" ../nvt.mdp ../npt.mdp ../md.mdp
grep "\[" index.ndx
```

b) Se `[ Protein_LMC ]` não existir no `index.ndx`, criar:

```bash
gmx make_ndx -f 4-solv_ions.gro -o index.ndx
> 1 | 13        # 13 = número do grupo LMC na lista
> q
```

c) Corrigir o `tc-grps` nos três `.mdp` (nvt, npt e md):

```bash
sed -i 's/^tc[-_]grps.*/tc-grps                 = Protein_LMC Water_and_ions/' ../nvt.mdp ../npt.mdp ../md.mdp
grep -n "tc-grps\|tau_t\|ref_t" ../nvt.mdp ../npt.mdp ../md.mdp
```

`tau_t` e `ref_t` precisam ter **2 valores** cada (um para cada grupo):

```
tc-grps   = Protein_LMC  Water_and_ions
tau_t     = 0.1          0.1
ref_t     = 300          300
```

d) Submeter de novo: `sbatch md_dk2.job`.

> O ligante é acoplado junto com a proteína porque, sozinho (poucos átomos), a temperatura dele ficaria instável.

### 5. Warning no `grompp` do NVT: 116 átomos fora dos grupos de remoção do centro de massa

```
WARNING 1 [file ../nvt.mdp]:
  116 atoms are not part of any center of mass motion removal group.
...
Fatal error:
Too many warnings (1).
```

**Causa:** mesmo problema do item 4, mas na linha `comm-grps` (remoção do movimento do centro de massa): ela lista grupos sem o ligante (ex.: `Protein Water_and_ions`).

**Solução:** usar um único grupo para o sistema todo (o que o próprio GROMACS recomenda). **Não** usar `-maxwarn`.

```bash
grep -n "comm[-_]grps" ../nvt.mdp ../npt.mdp ../md.mdp
sed -i 's/^comm[-_]grps.*/comm-grps               = System/' ../nvt.mdp ../npt.mdp ../md.mdp
grep -n "comm[-_]grps" ../nvt.mdp ../npt.mdp ../md.mdp
```

> A NOTE 2 (remoção do centro de massa com restrição de posição) é normal no NVT/NPT e pode ser ignorada — é só nota, não conta como warning. As mensagens `Ignoring obsolete mdp entry 'title'` / `'ns_type'` também são inofensivas.

### 6. Warning no `grompp` do NPT: barostato Berendsen

```
WARNING 1 [file ../npt.mdp]:
  The Berendsen barostat does not generate any strictly correct ensemble,
  ...we would recommend the C-rescale barostat...
Fatal error:
Too many warnings (1).
```

**Causa:** o `npt.mdp` usa `pcoupl = Berendsen`. No GROMACS 2022 isso gera warning e o `grompp` para.

**Solução:** trocar para `C-rescale` (disponível desde o GROMACS 2021). **Não** usar `-maxwarn`.

```bash
grep -n "pcoupl\|tau_p\|tau-p" ../npt.mdp ../md.mdp
sed -i 's/^\(pcoupl[ \t]*=[ \t]*\)[Bb]erendsen/\1C-rescale/' ../npt.mdp ../md.mdp
grep -n "pcoupl" ../npt.mdp ../md.mdp      # deve mostrar C-rescale
```

**Retomar do NPT sem refazer minimização e NVT** (se o NVT já terminou — `md_dk2.out` mostra "NVT OK" e existem `6-nvt.gro` e `6-nvt.cpt`):

```bash
cp md_dk2.job md_dk2_npt.job
nano md_dk2_npt.job
```

No `md_dk2_npt.job`: trocar `md_dk2` por `md_dk2_npt` nas linhas `--job-name`, `--output` e `--error`, e **apagar** os blocos "Minimização" e "Equilíbrio NVT" (deixar só NPT e Produção). Depois:

```bash
sbatch md_dk2_npt.job
```

### 7. Arquivos `#nome.1#` na pasta

São backups que o GROMACS cria quando um arquivo é sobrescrito (ex.: ao submeter o job de novo). Podem ser apagados:

```bash
rm \#*
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

### minim.mdp

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
