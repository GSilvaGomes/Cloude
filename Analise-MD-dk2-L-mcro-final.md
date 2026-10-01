# Análise da Dinâmica Molecular — FINAL (100 ns)

## POSE DOCKING - L_mcro_TRPA1_dk2

Sistema: TRPA1 + L-mcro (resíduo `LMC`), CHARMM36, TIP3P, NaCl 0,1 M, 100 ns de produção.
Pasta no Vital: `/storage/zuleika/volume2/project/gisele_picolo/minicro_docking/MD/MD_L_mcro_TRPA1_dk2`

> Roteiro completo, em ordem, para analisar a MD **depois que a produção terminar**.
> A análise parcial (~57 ns) está em `Analise-MD-dk2-L-mcro-parcial.md` — usar para comparar.
> Explicação de cada comando e opção: item 10 do arquivo parcial.

---

## Visão geral

| Passo | O que fazer | Onde |
|---|---|---|
| 0 | Conferir que a MD terminou e guardar os resultados parciais | Vital |
| 1 | Centralizar a trajetória | Vital |
| 2 | Conferir os grupos do índice e criar o grupo do bolsão | Vital |
| 3 | Definir a parte estável (platô) | Vital |
| 4 | RMSD: proteína, ligante, complexo, bolsão, ligante no bolsão, interno | Vital |
| 5 | RMSF: proteína e ligante | Vital |
| 6 | Distância ligante–bolsão | Vital |
| 7 | Ligações de hidrogênio e contatos | Vital |
| 8 | Médias ± desvio na parte estável | Vital |
| 9 | Frames e pose representativa (cluster) | Vital |
| 10 | Copiar os resultados para o Ubuntu | Ubuntu |
| 11 | Gerar os gráficos (gnuplot → PDF) | Ubuntu |
| 12 | Preencher a tabela de resultados e comparar com o parcial | — |

Todos os comandos dos passos 0 a 9 são rodados **dentro da pasta do sistema no Vital**.

---

## Passo 0 — Conferir que a MD terminou e guardar o parcial

```bash
cat md_dk2.out                 # deve ter "Produção OK"
tail -3 8-md.log               # deve ter "Finished mdrun"
grep "Performance" 8-md.log    # velocidade final (ns/dia)
```

Guardar os resultados parciais (os arquivos `9-*` serão recriados):

```bash
mkdir -p parcial_57ns
mv 9-* frame_*ps.pdb parcial_57ns/ 2>/dev/null
ls parcial_57ns/
```

---

## Passo 1 — Centralizar a trajetória

```bash
echo "Protein_LMC System" | gmx trjconv -s 8-md.tpr -f 8-md.xtc -n index.ndx -o 9-md_center.xtc -center -pbc mol -ur compact
```

Trajetória leve para visualizar (só proteína + ligante, 1 a cada 10 frames) e estrutura de referência:

```bash
echo "24 24" | gmx trjconv -s 8-md.tpr -f 9-md_center.xtc -n index.ndx -o 9-md_vis.xtc -skip 10
echo "24 24" | gmx trjconv -s 8-md.tpr -f 9-md_center.xtc -n index.ndx -o 9-ref.pdb -dump 0
```

---

## Passo 2 — Grupos do índice e grupo do bolsão

O nome `LMC` aparece duas vezes no `index.ndx` → usar sempre **números**. Conferir:

```bash
grep "\[" index.ndx | awk '{print NR-1, $2}'
```

| Grupo | Número |
|---|---|
| Protein | **1** |
| C-alpha | **3** |
| Backbone | **4** |
| LMC (ligante, 116 átomos) | **13** |
| Protein_LMC | **24** |

Criar o grupo do bolsão (backbone dos resíduos a até 5 Å do ligante no início da produção):

```bash
rm -f pocket.ndx index_pocket.ndx
gmx select -s 8-md.tpr -n index.ndx -select 'group "Backbone" and same residue as within 0.5 of resname LMC' -on pocket.ndx
cat index.ndx pocket.ndx > index_pocket.ndx
grep "\[" index_pocket.ndx | awk '{print NR-1, $2}' | tail -2
```

O último grupo da lista é o bolsão (**25**). Se for outro número, trocar `25` nos passos seguintes.

---

## Passo 3 — Definir a parte estável (platô)

Primeiro calcular o RMSD da proteína (passo 4.1), olhar o gráfico e escolher a partir de quantos ns a curva estabiliza. Na análise parcial, o ligante estabilizou depois de ~25 ns.

Definir as variáveis (trocar `30` pelo início do platô escolhido):

```bash
T0=30                 # início do platô em ns (para as médias)
B=$((T0*1000))        # o mesmo em ps (para os comandos com -b)
echo "Platô a partir de $T0 ns ($B ps)"
```

> As variáveis valem só nesta sessão do terminal. Se fechar o terminal, definir de novo.

---

## Passo 4 — RMSD

```bash
# 4.1 Proteína (backbone)
echo "4 4"   | gmx rms -s 8-md.tpr -f 9-md_center.xtc -n index.ndx        -o 9-rmsd_prot.xvg        -tu ns
# 4.2 Ligante (ajuste no backbone inteiro)
echo "4 13"  | gmx rms -s 8-md.tpr -f 9-md_center.xtc -n index.ndx        -o 9-rmsd_lig.xvg         -tu ns
# 4.3 Complexo
echo "4 24"  | gmx rms -s 8-md.tpr -f 9-md_center.xtc -n index.ndx        -o 9-rmsd_complexo.xvg    -tu ns
# 4.4 Ligante (ajuste no bolsão) — o mais importante para o ligante
echo "25 13" | gmx rms -s 8-md.tpr -f 9-md_center.xtc -n index_pocket.ndx -o 9-rmsd_lig_pocket.xvg  -tu ns
# 4.5 Bolsão
echo "25 25" | gmx rms -s 8-md.tpr -f 9-md_center.xtc -n index_pocket.ndx -o 9-rmsd_pocket.xvg      -tu ns
# 4.6 Interno do ligante (conformação)
echo "13 13" | gmx rms -s 8-md.tpr -f 9-md_center.xtc -n index.ndx        -o 9-rmsd_lig_interno.xvg -tu ns

# Converter tudo para Å
for f in prot lig complexo lig_pocket pocket lig_interno; do awk '!/^[#@]/{print $1, $2*10}' 9-rmsd_$f.xvg > 9-rmsd_${f}_A.dat; done
```

---

## Passo 5 — RMSF

```bash
# 5.1 Proteína, por resíduo (C-alpha), na parte estável
echo "3" | gmx rmsf -s 8-md.tpr -f 9-md_center.xtc -n index.ndx -o 9-rmsf.xvg -res -b $B -oq 9-rmsf_prot_bfac.pdb
awk '!/^[#@]/{print $1, $2*10}' 9-rmsf.xvg > 9-rmsf_A.dat

# 5.2 Ligante — flexibilidade interna (por átomo)
echo "13" | gmx rmsf -s 8-md.tpr -f 9-md_center.xtc -n index.ndx -o 9-rmsf_lig.xvg -b $B -oq 9-rmsf_lig_bfac.pdb
awk '!/^[#@]/{print $1, $2*10}' 9-rmsf_lig.xvg > 9-rmsf_lig_A.dat

# 5.3 Ligante — movimento dentro do sítio (trajetória alinhada pelo bolsão)
echo "25 0" | gmx trjconv -s 8-md.tpr -f 9-md_center.xtc -n index_pocket.ndx -o 9-md_fit_pocket.xtc -fit rot+trans
echo "13" | gmx rmsf -s 8-md.tpr -f 9-md_fit_pocket.xtc -n index.ndx -o 9-rmsf_lig_pocket.xvg -nofit -b $B -oq 9-rmsf_lig_pocket_bfac.pdb
awk '!/^[#@]/{print $1, $2*10}' 9-rmsf_lig_pocket.xvg > 9-rmsf_lig_pocket_A.dat
```

Os `*_bfac.pdb` têm o RMSF na coluna B-factor → ver em cores (VMD: *Coloring Method → Beta*; PyMOL: `spectrum b, blue_red`).

---

## Passo 6 — Distância ligante–bolsão

```bash
gmx distance -s 8-md.tpr -f 9-md_center.xtc -n index_pocket.ndx -select 'com of group 13 plus com of group 25' -oall 9-dist_lig_pocket.xvg -tu ns
awk '!/^[#@]/{print $1, $2*10}' 9-dist_lig_pocket.xvg > 9-dist_lig_pocket_A.dat
```

---

## Passo 7 — Ligações de hidrogênio e contatos

```bash
# 7.1 Ligações de H proteína–ligante
echo "1 13" | gmx hbond -s 8-md.tpr -f 9-md_center.xtc -n index.ndx -num 9-hbond.xvg -tu ns
awk '!/^[#@]/{print $1, $2}' 9-hbond.xvg > 9-hbond.dat

# 7.2 Contatos (< 4 Å) e distância mínima
echo "1 13" | gmx mindist -s 8-md.tpr -f 9-md_center.xtc -n index.ndx -od 9-mindist.xvg -on 9-contatos.xvg -d 0.4 -tu ns
awk '!/^[#@]/{print $1, $2}'    9-contatos.xvg > 9-contatos.dat
awk '!/^[#@]/{print $1, $2*10}' 9-mindist.xvg  > 9-mindist_A.dat
```

---

## Passo 8 — Médias ± desvio (parte estável)

```bash
echo "=== Médias a partir de $T0 ns ===" | tee 9-medias.txt
for f in prot lig complexo lig_pocket pocket lig_interno; do awk -v t0=$T0 -v n=$f '$1>=t0{s+=$2; q+=$2*$2; c++} END{m=s/c; printf "RMSD %-12s = %5.2f ± %4.2f Å\n", n, m, sqrt(q/c-m*m)}' 9-rmsd_${f}_A.dat; done | tee -a 9-medias.txt
awk -v t0=$T0 '$1>=t0{s+=$2; q+=$2*$2; c++} END{m=s/c; printf "Distância lig-bolsão = %5.2f ± %4.2f Å\n", m, sqrt(q/c-m*m)}' 9-dist_lig_pocket_A.dat | tee -a 9-medias.txt
for f in hbond contatos; do awk -v t0=$T0 -v n=$f '$1>=t0{s+=$2; q+=$2*$2; c++} END{m=s/c; printf "%-20s = %6.1f ± %5.1f\n", n, m, sqrt(q/c-m*m)}' 9-$f.dat; done | tee -a 9-medias.txt
awk -v t0=$T0 '$1>=t0{s+=$2; q+=$2*$2; c++} END{m=s/c; printf "Distância mínima     = %5.2f ± %4.2f Å\n", m, sqrt(q/c-m*m)}' 9-mindist_A.dat | tee -a 9-medias.txt
```

Os valores ficam salvos em **`9-medias.txt`**.

---

## Passo 9 — Frames e pose representativa

### 9.1 — Frames em tempos-chave

Tempos em ps. Os da análise parcial foram ~4, ~19 e ~38 ns; acrescentar outros tempos se aparecerem novos saltos nos gráficos.

```bash
for t in 0 4000 19000 38000 50000 75000 100000; do echo "24" | gmx trjconv -s 8-md.tpr -f 9-md_center.xtc -n index.ndx -o frame_${t}ps.pdb -dump $t; done
```

### 9.2 — Pose representativa do ligante (análise de clusters) — opcional

Agrupa as poses do ligante na parte estável (trajetória alinhada pelo bolsão do passo 5.3). Corte de 2 Å (0,2 nm).

```bash
echo "13 24" | gmx cluster -s 8-md.tpr -f 9-md_fit_pocket.xtc -n index.ndx -method gromos -cutoff 0.2 -nofit -b $B -skip 5 -g 9-cluster.log -cl 9-clusters.pdb -sz 9-cluster_size.xvg
grep -A5 "cl. | #st" 9-cluster.log
```

- Grupo 13 = LMC (para comparar as poses); grupo 24 = Protein_LMC (o que sai no `.pdb`).
- `9-clusters.pdb`: a estrutura central de cada cluster. O **cluster 1** é a pose mais frequente → usar para figuras.
- `-skip 5`: usa 1 a cada 5 frames (o cálculo cresce muito com o número de frames).

### 9.3 — Resíduos do sítio na pose representativa (PyMOL)

```
load 9-clusters.pdb
split_states 9-clusters
select sitio, byres (polymer.protein within 4 of resn LMC) and 9-clusters_0001
show sticks, sitio
iterate sitio and name CA, print(resn, resi)
```

---

## Passo 10 — Copiar para o Ubuntu

No **terminal do Ubuntu**, numa pasta nova para os resultados finais:

```bash
mkdir -p ~/MD_dk2_final && cd ~/MD_dk2_final
V="geniana_gomes@endereco_do_vital:/storage/zuleika/volume2/project/gisele_picolo/minicro_docking/MD/MD_L_mcro_TRPA1_dk2"
scp "$V/9-*_A.dat" "$V/9-hbond.dat" "$V/9-contatos.dat" "$V/9-medias.txt" .
scp "$V/9-*_bfac.pdb" "$V/frame_*ps.pdb" "$V/9-clusters.pdb" "$V/9-ref.pdb" "$V/9-md_vis.xtc" .
ls
```

Trocar `endereco_do_vital` pelo endereço usado no `ssh`.

---

## Passo 11 — Gráficos (gnuplot → PDF)

Criar cada script com `nano nome.gp`, colar, salvar (`Ctrl + O`, `Enter`, `Ctrl + X`) e rodar `gnuplot nome.gp`.

### 11.1 — `plot_rmsd.gp` → `rmsd_final_dk2.pdf`

```gnuplot
set terminal pdfcairo size 10,5 font "Arial,12" enhanced
set output "rmsd_final_dk2.pdf"
set title "RMSD - TRPA1 + L-mcro (dk2) - 100 ns"
set xlabel "Tempo (ns)"
set ylabel "RMSD (Å)"
set xrange [0:100]
set grid
set key top left
plot "9-rmsd_prot_A.dat"     with lines lw 2 lc rgb "red"   title "Proteína (backbone)", \
     "9-rmsd_lig_A.dat"      with lines lw 2 lc rgb "blue"  title "Ligante (LMC)", \
     "9-rmsd_complexo_A.dat" with lines lw 2 lc rgb "black" title "Complexo"
```

### 11.2 — `plot_rmsf.gp` → `rmsf_final_dk2.pdf`

```gnuplot
set terminal pdfcairo size 10,5 font "Arial,12" enhanced
set output "rmsf_final_dk2.pdf"
set title "RMSF por resíduo (C-alpha) - TRPA1 + L-mcro (dk2) - 100 ns"
set xlabel "Resíduo"
set ylabel "RMSF (Å)"
set grid
unset key
plot "9-rmsf_A.dat" with lines lw 2 lc rgb "red"
```

### 11.3 — `plot_lig_sitio.gp` → `ligante_sitio_final_dk2.pdf`

```gnuplot
set terminal pdfcairo size 10,8 font "Arial,12" enhanced
set output "ligante_sitio_final_dk2.pdf"
set multiplot layout 2,1 title "L-mcro no sítio da TRPA1 (dk2) - 100 ns"

set xlabel "Tempo (ns)"
set xrange [0:100]
set grid
set ylabel "RMSD (Å)"
set key top left
plot "9-rmsd_lig_A.dat"         with lines lw 1.5 lc rgb "#9ecae1" title "Ligante (ajuste: backbone inteiro)", \
     "9-rmsd_lig_pocket_A.dat"  with lines lw 2   lc rgb "blue"    title "Ligante (ajuste: bolsão)", \
     "9-rmsd_pocket_A.dat"      with lines lw 2   lc rgb "red"     title "Bolsão (backbone)", \
     "9-rmsd_lig_interno_A.dat" with lines lw 2   lc rgb "#2ca02c" title "Ligante (conformação interna)"

set ylabel "Distância (Å)"
unset key
plot "9-dist_lig_pocket_A.dat" with lines lw 2 lc rgb "black"

unset multiplot
```

### 11.4 — `plot_rmsf_lig.gp` → `rmsf_ligante_final_dk2.pdf`

```gnuplot
set terminal pdfcairo size 10,5 font "Arial,12" enhanced
set output "rmsf_ligante_final_dk2.pdf"
set title "RMSF do ligante por átomo - L-mcro (dk2) - 100 ns"
set xlabel "Átomo (número no sistema)"
set ylabel "RMSF (Å)"
set grid
set key top left
plot "9-rmsf_lig_A.dat"        with linespoints lw 2 pt 7 ps 0.5 lc rgb "blue"    title "Flexibilidade interna", \
     "9-rmsf_lig_pocket_A.dat" with linespoints lw 2 pt 7 ps 0.5 lc rgb "#9ecae1" title "Movimento no sítio (ajuste: bolsão)"
```

### 11.5 — `plot_interacoes.gp` → `interacoes_final_dk2.pdf`

```gnuplot
set terminal pdfcairo size 10,8 font "Arial,12" enhanced
set output "interacoes_final_dk2.pdf"
set multiplot layout 2,1 title "Interações L-mcro – TRPA1 (dk2) - 100 ns"

set xlabel "Tempo (ns)"
set xrange [0:100]
set grid
unset key
set ylabel "Ligações de H"
plot "9-hbond.dat" with lines lw 1.5 lc rgb "blue"

set ylabel "Contatos (< 4 Å)"
plot "9-contatos.dat" with lines lw 1.5 lc rgb "red"

unset multiplot
```

### 11.6 — Gerar todos

```bash
for g in plot_rmsd plot_rmsf plot_lig_sitio plot_rmsf_lig plot_interacoes; do gnuplot $g.gp && echo "$g OK"; done
ls *_final_dk2.pdf
```

Abrir: dois cliques no gerenciador de arquivos, ou `explorer.exe arquivo.pdf` (WSL).

### 11.7 — Visualizar a trajetória final (VMD)

```bash
vmd 9-ref.pdb 9-md_vis.xtc
```

No **Tk Console**: proteína em cartoon vermelho, ligante em sticks azul sem H:

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

---

## Passo 12 — Resultados e comparação com o parcial

### 12.1 — Tabela (média ± desvio na parte estável)

Copiar os valores do `9-medias.txt`.

| Medida | Parcial (25–57 ns) | Final (___–100 ns) |
|---|---|---|
| RMSD proteína (backbone) | 12,01 ± 1,05 Å | |
| RMSD complexo | — | |
| RMSD ligante (ajuste: backbone inteiro) | 7,41 ± 0,78 Å | |
| RMSD ligante (ajuste: bolsão) | 5,71 ± 0,81 Å | |
| RMSD bolsão (backbone) | 1,86 ± 0,43 Å | |
| RMSD interno do ligante | 3,07 ± 0,33 Å | |
| Distância ligante–bolsão | 6,23 ± 0,57 Å | |
| Ligações de H (nº) | — | |
| Contatos < 4 Å (nº) | — | |
| Distância mínima | — | |
| Velocidade da produção (ns/dia) | — | |

### 12.2 — Perguntas para responder com os dados finais

1. O RMSD da proteína **estabilizou** (platô) até 100 ns? Em quanto?
2. O **bolsão** continuou estável (~1–2 Å)?
3. O ligante **permaneceu no sítio** até o fim (distância estável, sem tendência de aumento)?
4. A **pose** encontrada depois de ~20 ns (RMSD no bolsão ~5,7 Å) se manteve, ou houve nova mudança?
5. Quantas **ligações de H** e **contatos** em média? Algum período sem interação?
6. Quais **resíduos** formam o sítio na pose representativa (cluster 1)?
7. O cluster 1 representa que **fração** das poses (estabilidade da pose)?

### 12.3 — Valores de referência para interpretar

| Medida | Estável | Atenção |
|---|---|---|
| RMSD proteína | Platô (~1–3 Å; 3–4 Å em proteínas grandes/flexíveis) | Subindo sem estabilizar |
| RMSD bolsão | ~1–2 Å | > 3 Å → sítio deformado |
| RMSD ligante (ajuste: bolsão) | < 2–3 Å (mesma pose do docking) | > 4–5 Å ou saltos → mudou de pose |
| Distância ligante–bolsão | Constante | Aumentando → saindo do sítio |
| Ligações de H / contatos | Estáveis, sem cair a zero | Quedas longas → perdendo interação |
| RMSF | Picos em terminais/alças | Picos nos resíduos do sítio |

> Lembrar das particularidades do sistema (ver item 11.1 do parcial): canal de membrana simulado só em água, uma cadeia do tetrâmero, modelo AlphaFold → RMSD global da proteína alto é esperado; o mais informativo para o ligante são o **bolsão** e o **ligante com ajuste no bolsão**.

### 12.4 — Conclusão final

_(escrever depois de preencher a tabela)_
