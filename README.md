# hse2025

Cimparative analysis of Lactobacillus helveticus strains derived from milk- and cheese-associated samples


### Task 1. What Makes the Tajik Strain Unique?

```
import pandas as pd

# Load Roary gene presence matrix
df_roary = pd.read_csv("gene_presence_absence.csv")

# Tajik strain
tajik_strain = "Lactobacillus_helveticus_TJA10"
other_strains = [
    'Lactobacillus_helveticus_DPC_4571',
    'Lactobacillus_helveticus_H10',
    'Lactobacillus_helveticus_CNRZ32',
    'Lactobacillus_helveticus_Lh21462',
    'Lactobacillus_helveticus_Lh21463',
    'Lactobacillus_helveticus_Lh11961',
    'Lactobacillus_helveticus_Lh11051',
    'Lactobacillus_helveticus_Lh21456'
]

# Filter for unique Tajik genes
tajik_unique = df_roary[
    (df_roary[tajik_strain] == tajik_strain) & 
    (df_roary[other_strains].isna().all(axis=1))
]

# Shared among all others, but missing in Tajik strain
shared_others = df_roary[
    (df_roary[tajik_strain].isna()) &
    (df_roary[other_strains].notna().all(axis=1))
]

# Output summary
print("Unique genes in Tajik strain only:", tajik_unique.shape[0])
print("Genes shared by others but missing in Tajik strain:", shared_others.shape[0])

# Save to CSV
tajik_unique.to_csv("tajik_unique_genes.csv", index=False)
shared_others.to_csv("tajik_missing_core.csv", index=False)
```

Gene presence/absence heatmap:
```
import seaborn as sns
import matplotlib.pyplot as plt

heatmap_df = df_roary[[tajik_strain] + other_strains].notna().astype(int).T
sns.heatmap(heatmap_df, cbar=False)
plt.title("Gene Presence/Absence Heatmap")
plt.xlabel("Gene Clusters")
plt.ylabel("Strains")
plt.show()
```

Get protein sequences from prokka fasta file
```
from Bio import SeqIO

# Load Tajik strain Prokka proteins
faa_file = "prokka_outputs_helvet/Lactobacillus_helveticus_TJA10/Lactobacillus_helveticus_TJA10.faa"
tajik_gene_ids = ["PROKKA_00001", "PROKKA_00034", ...]  # fill from Roary

records = SeqIO.parse(faa_file, "fasta")
filtered = [r for r in records if any(gid in r.id for gid in tajik_gene_ids)]

# Save for EggNOG
SeqIO.write(filtered, "tajik_specific_proteins.faa", "fasta")

```



