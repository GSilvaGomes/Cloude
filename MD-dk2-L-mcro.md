# Dinâmica Molecular

## POSE DOCKING - L_mcro_TRPA1_dk2

## Script construção arquivos para MD

- Campo de força: CHARMM36 (2020)
- Caixa de água: 1.2 nm, modelo TIP3P
- Sais: NaCl [0,1 M]

## Preparo topologia ligante

Topologia preparada pela plataforma CGenFF (https://cgenff.com/) para o campo de força CHARMM36 e convertida para o GROMACS.

Arquivos gerados:

```
L-mcro_dk2_gmx.pdb
L-mcro_dk2_gmx.itp
```

## Protonação pós-docking

Programa PDB2PQR (https://server.poissonboltzmann.org/pdb2pqr), que gera os arquivos `.log` e `.pqr`.

Script Python para saber quais aminoácidos da proteína devem estar protonados (pH 7):

```bash
python3 Output-ProtonationPDB2PQR.py TRPA1_AF3_Y_dk2.log TRPA1_AF3_Y_dk2.pqr 7 > output_TRPA1_AF3_Y_dk2.txt
```

## Topologia proteína

```bash
gmx pdb2gmx -f TRPA1_AF3_Y_dk2.pdb -water tip3p -ignh -lys -his -asp -glu -ter -o conf.pdb
```

