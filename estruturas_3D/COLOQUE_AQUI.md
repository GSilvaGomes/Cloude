# Estruturas 3D para conferencia

Coloque nesta pasta os arquivos 3D dos ligantes, um por composto, com o nome
do composto no arquivo: `mcro.sdf`, `mcr1.sdf`, ... ou `.pdb` / `.mol2`.

Depois rode, da raiz do projeto:

```bash
python3 scripts/checa_quiralidade_3D.py estruturas_3D/
```

O script confere, para cada estrutura:

- quantos centros quirais tem (esperado: 6) e o CIP de cada um;
- se todos sao S, ou seja, all-L -- um centro R significa residuo D;
- se a formula molecular bate com o SMILES de referencia;
- se a conectividade bate com o SMILES de referencia.

## Formato

Prefira **SDF** ou **MOL2**: guardam ordem de ligacao explicita.

**PDB funciona, mas e o menos confiavel** -- nao guarda ordem de ligacao,
entao a percepcao pode falhar em residuo fora do padrao, como o pGlu. Se o
PDB der erro de leitura, exporte SDF do mesmo arquivo e rode de novo.
