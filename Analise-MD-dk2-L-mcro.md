# Análise da Dinâmica Molecular

## POSE DOCKING - L_mcro_TRPA1_dk2

Sistema: TRPA1 + L-mcro (resíduo `LMC`), CHARMM36, TIP3P, NaCl 0,1 M.
Pasta no Vital: `/storage/zuleika/volume2/project/gisele_picolo/minicro_docking/MD/MD_L_mcro_TRPA1_dk2`

Todos os comandos são rodados **dentro da pasta do sistema**. Nenhum deles para a MD — só leem os arquivos.

---

## 1. Acompanhar a MD enquanto roda

**Etapa atual e previsão de término:**

```bash
for f in 5-minim 6-nvt 7-npt 8-md; do [ -f ${f}_run.log ] && echo "$f: $(tr '\r' '\n' < ${f}_run.log | grep -E 'will finish|Step=|^step' | tail -1)"; done
```

**Quantos ns da produção já foram simulados** (segundo número = tempo em ps):

```bash
grep -A1 "Step           Time" 8-md.log | tail -1
```

**Job ainda rodando? Quais etapas já terminaram?**

```bash
squeue -u $USER      # ST: R = rodando | CD = terminou | F = falhou
cat md_dk2.out       # Minimização OK, NVT OK, NPT OK, Produção OK
```

**Atualizando sozinho a cada 60 s** (`Ctrl + C` fecha só a visualização):

```bash
watch -n 60 "tr '\r' '\n' < 8-md_run.log | grep 'will finish' | tail -1"
```

**Velocidade da simulação** (no fim de cada etapa concluída; 1º número = ns/dia):

```bash
grep -H "Performance" 6-nvt.log 7-npt.log 8-md.log
```

> Referência: NVT anterior rodou a ~80 ns/dia (0,30 h/ns) com 1 GPU → 100 ns de produção ≈ 30 h.

---

## 2. Visualizar a trajetória (VMD / PyMOL)

### 2.1 — No Vital: trajetória centralizada, só proteína + ligante

Funciona com a MD ainda rodando (usa os frames que já existem).

```bash
echo "Protein_LMC Protein_LMC" | gmx trjconv -s 8-md.tpr -f 8-md.xtc -n index.ndx -o parcial.xtc -center -pbc mol -ur compact -skip 10
echo "Protein_LMC Protein_LMC" | gmx trjconv -s 8-md.tpr -f 8-md.xtc -n index.ndx -o ref.pdb -center -pbc mol -ur compact -dump 0
```

- `-skip 10`: 1 a cada 10 frames (arquivo mais leve). Tirar para ver todos.
- Para atualizar depois, rodar de novo só o primeiro comando.

### 2.2 — Baixar para o computador

No **terminal do seu computador**:

```bash
scp geniana_gomes@endereco_do_vital:/storage/zuleika/volume2/project/gisele_picolo/minicro_docking/MD/MD_L_mcro_TRPA1_dk2/{ref.pdb,parcial.xtc} .
```

### 2.3 — VMD

```bash
vmd ref.pdb parcial.xtc
```

Em **Extensions → Tk Console**, colar:

**Proteína em cartoon vermelho + ligante em sticks azul, sem hidrogênios:**

```tcl
while {[molinfo top get numreps] > 0} { mol delrep 0 top }

mol representation NewCartoon
mol selection "protein"
mol color ColorID 1
mol addrep top

mol representation Licorice 0.3 12 12
mol selection "resname LMC and not hydrogen"
mol color ColorID 0
mol addrep top
```

**Opcional — resíduos da proteína a até 5 Å do ligante** (atualiza a cada frame):

```tcl
mol representation Licorice 0.2 12 12
mol selection "protein and same residue as within 5 of resname LMC"
mol color Name
mol addrep top
mol selupdate [expr [molinfo top get numreps] - 1] top 1
```

Cores (`ColorID`): 0 azul · 1 vermelho · 4 amarelo · 7 verde · 10 ciano · 21 azul-escuro · 29 vermelho-escuro.

Pelo menu: **Graphics → Representations** → *Selected Atoms*, *Drawing Method*, *Coloring Method* → *Apply* (*Create Rep* para uma nova).

### 2.4 — PyMOL

```
load ref.pdb, complexo
load_traj parcial.xtc, complexo
hide everything
show cartoon, polymer.protein
color red, polymer.protein
show sticks, resn LMC and not elem H
color blue, resn LMC and elem C
```

Usar os botões de play no canto inferior direito.

---

## 3. RMSD da proteína, do ligante e do complexo (em Å)

### 3.1 — Centralizar a trajetória (uma vez)

```bash
echo "Protein_LMC System" | gmx trjconv -s 8-md.tpr -f 8-md.xtc -n index.ndx -o 9-md_center.xtc -center -pbc mol -ur compact
```

### 3.2 — Calcular o RMSD

Ajuste (*fit*) sempre pelo **backbone da proteína** → o RMSD do ligante mostra quanto ele se moveu **dentro do sítio**.

```bash
# Proteína (backbone)
echo "Backbone Backbone" | gmx rms -s 8-md.tpr -f 9-md_center.xtc -n index.ndx -o 9-rmsd_prot.xvg -tu ns

# Ligante (LMC), em relação à proteína
echo "Backbone LMC" | gmx rms -s 8-md.tpr -f 9-md_center.xtc -n index.ndx -o 9-rmsd_lig.xvg -tu ns

# Complexo (proteína + ligante)
echo "Backbone Protein_LMC" | gmx rms -s 8-md.tpr -f 9-md_center.xtc -n index.ndx -o 9-rmsd_complexo.xvg -tu ns
```

### 3.3 — Converter nm → Å (× 10)

```bash
for f in prot lig complexo; do awk '!/^[#@]/{print $1, $2*10}' 9-rmsd_$f.xvg > 9-rmsd_${f}_A.dat; done
```

Gera `9-rmsd_prot_A.dat`, `9-rmsd_lig_A.dat`, `9-rmsd_complexo_A.dat` (colunas: tempo em ns | RMSD em Å).

### 3.4 — Média e desvio padrão (parte estável)

Trocar `20` pelo tempo (ns) a partir do qual a curva estabiliza:

```bash
for f in prot lig complexo; do awk -v t0=20 -v n=$f '$1>=t0{s+=$2; q+=$2*$2; c++} END{m=s/c; printf "%-9s média = %.2f Å  ±  %.2f Å  (%d frames)\n", n, m, sqrt(q/c-m*m), c}' 9-rmsd_${f}_A.dat; done
```

---

## 4. Gráficos com gnuplot

### 4.1 — Instalar

- **Ubuntu (seu computador):**
  ```bash
  sudo apt update
  sudo apt install gnuplot
  ```
- **Vital:** verificar com `which gnuplot`. Se não tiver (e sem `sudo`), baixar os `.dat` com `scp` e plotar no seu computador.
- **Windows:** https://sourceforge.net/projects/gnuplot/files/gnuplot/

Alternativa: `sudo apt install grace` → `xmgrace 9-rmsd_prot.xvg 9-rmsd_lig.xvg` (abre o `.xvg` direto, mas em nm).

### 4.2 — Script do gráfico

```bash
nano plot_rmsd.gp
```

Colar:

```gnuplot
set terminal pngcairo size 1400,700 font "Arial,16"
set output "rmsd_dk2.png"
set xlabel "Tempo (ns)"
set ylabel "RMSD (Å)"
set title "RMSD - TRPA1 + L-mcro (dk2)"
set key top left
set grid
plot "9-rmsd_prot_A.dat"     with lines lw 2 lc rgb "red"   title "Proteína (backbone)", \
     "9-rmsd_lig_A.dat"      with lines lw 2 lc rgb "blue"  title "Ligante (LMC)", \
     "9-rmsd_complexo_A.dat" with lines lw 2 lc rgb "black" title "Complexo"
```

Salvar (`Ctrl + O`, `Enter`, `Ctrl + X`) e gerar a imagem:

```bash
gnuplot plot_rmsd.gp      # gera rmsd_dk2.png
```

Para tirar uma curva, apagar a linha dela no `plot` (e a vírgula `, \` da linha anterior).

---

## 5. Como interpretar o RMSD

| Curva | Esperado | Alerta |
|---|---|---|
| **Proteína** | Sobe nos primeiros ns e forma **platô** (~1–3 Å) | Sobe sem parar → estrutura ainda mudando (comum em modelo AlphaFold) |
| **Ligante** | Baixo e estável (~1–3 Å) → **permaneceu no sítio** do docking | Saltos grandes (> 4–5 Å) → mudou de pose ou saiu do sítio |
| **Complexo** | Próximo da proteína (ela tem muito mais átomos) | — |

- Usar a **parte estável (platô)** para médias e análises; os primeiros ns ainda são equilíbrio.
- Se o ligante der um salto, abrir a trajetória no VMD naquele tempo para ver o que aconteceu.
