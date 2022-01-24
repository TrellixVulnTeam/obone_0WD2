
# Bone Package

Python workflow for boolean lab analysis.


## Features

- NCBI GEO parser to build expression and survival files
- Preprocessing to read raw tarfile



## Demo

```
import matplotlib.pyplot as plt

import bone

# create the survival and expression file with the GEO class
gse40060 = bone.GEO(accessionID="GSE40060")
gse40060_survival = gse40060.survival()  # if unspecified, first GPL is used
gse40060_expr = gse40060.expr(rename_genes=True, probeID="ENSG")
gse40060_expr = gse40060_expr.fillna(0)

# initialize bone with expression and survival file
my_bone = bone.BoNE(expr=gse40060_expr, survival=gse40060_survival)

# define variables for bone
survival_name = "c source_name_ch1"
gse40060_groups = {"E": "endogenous", "O": "overexpressed"}
alz_gene_weights = {
    -3: ["SVOP", "CACNG3", "PCYOX1L", "BEX1", "TUBB3", "NRN1", "GAP43", "RGS4"],
    1: ["BGN", "EHD2", "FCGRT", "NT5DC2", "ITGB5", "PDGFRB", "GPR4", "LAMB2"],
}
my_bone.init(
    survival_col=survival_name, 
    gene_weights=alz_gene_weights, 
    groups=gse40060_groups,
)

# visualize violin plot
plt.figure(figsize=(10, 5), dpi=100)
my_bone.violin()
plt.show()
```

## Authors

- [@otucher](https://www.github.com/otucher)


