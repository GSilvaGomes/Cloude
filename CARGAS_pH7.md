# Estado de protonação e carga líquida em pH 7,4

Referência para padronizar a preparação dos 24 ligantes antes do docking.

## Grupos ionizáveis da série

Só a posição 5 é ionizável. Todo o resto é neutro em pH 7,4:

| grupo | posição | estado em pH 7,4 | carga |
|---|---|---|---|
| pGlu (lactama) | 1 | neutro — o N está na lactama, não é amina livre | 0 |
| N-terminal | — | bloqueado pelo pGlu, não existe como amina livre | 0 |
| Phe / Val / Leu | 2 | neutro | 0 |
| Ser / Pro | 3 | neutro | 0 |
| Pro / Leu | 4 | neutro | 0 |
| **Glu** | **5** | **carboxilato, COO⁻** | **−1** |
| **Lys** | **5** | **amônio, NH₃⁺** | **+1** |
| Asn | 6 | neutro (amida da cadeia lateral) | 0 |
| C-terminal | — | amidado, não é ácido livre | 0 |

## Consequência: a série se divide em dois grupos de carga

| carga | n | compostos |
|---|---|---|
| **−1** (Glu na pos 5) | 12 | mcro, mcr1, mcr3, mcr5, mcr8, mcr9, mcr10, mcr11, mcr14, mcr15, mcr18, mcr19 |
| **+1** (Lys na pos 5) | 12 | mcr2, mcr4, mcr6, mcr7, mcr12, mcr13, mcr16, mcr17, mcr20, mcr21, mcr22, mcr23 |

A troca Glu→Lys na posição 5 muda a carga líquida em **+2 unidades**.

## O que isso implica para o docking

Comparar score entre o subgrupo −1 e o subgrupo +1 é comparar moléculas de
carga líquida diferente. Se o sítio do TRPA1 for eletrostaticamente
polarizado, o termo coulômbico tende a dominar a diferença, e o ranking pode
refletir a carga em vez da complementaridade de forma que se quer medir.

Isso não é erro de construção — é uma propriedade real dos compostos, e
provavelmente parte do que se quer testar. Mas convém:

1. **Comparar dentro de cada subgrupo** (os 12 com Glu entre si, os 12 com
   Lys entre si), onde a carga é constante e as diferenças vêm das posições
   2, 3 e 4;
2. Ao comparar entre subgrupos, tratar a diferença de carga explicitamente,
   e não como se fosse só mais uma substituição;
3. **Padronizar a protonação de todos os 24** na preparação — mesmo pH, mesma
   ferramenta. Misturar Glu como COOH neutro em uns e COO⁻ em outros
   introduz uma diferença de carga que não é química, é de preparo.

## Verificação

Confira a carga de cada estrutura preparada contra a tabela acima. No PDB, a
carga formal fica nas colunas 79-80.

```bash
grep -E '^(ATOM|HETATM)' arquivo.pdb | cut -c79-80 | grep -v '^ *$' | sort | uniq -c
```
