
# Bone Package

Python workflow for boolean lab analysis.


## Features

- NCBI GEO parser to build expression and survival files
- Preprocessing to read raw tarfile



## Demo

```
import bone

gse40060 = bone.GEO(accessionID="GSE40060")
# use default GPL number for survival and expr
survival = gse40060.survival()
expr = gse40060.expr(rename_genes=True)
expr = expr.fillna(0)
expr = bone.add_probeID(expr, "Homo Sapiens", "ENSG")
my_bone = bone.BoNE(expr, survival)

# define variables for bone
name = "c source_name_ch1"
groups = {"E": "endogenous", "O": "overexpressed"}
gene_weights = {-3: ['SVOP', 'CACNG3', 'PCYOX1L', 'BEX1']
                 1: ['BGN', 'EHD2', 'FCGRT', 'NT5DC2', 'ITGB5']}                 
my_bone.init(survival_col = name, 
             gene_weights = gene_weights, 
             groups = groups)

# visualize violin plot
plt.figure(figsize=(10, 5), dpi=100)
my_bone.violin()
plt.show()
```

## Authors

- [@otucher](https://www.github.com/otucher)


