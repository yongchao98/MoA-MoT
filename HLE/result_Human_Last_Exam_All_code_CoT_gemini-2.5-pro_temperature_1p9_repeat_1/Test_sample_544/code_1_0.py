import sys
# This script determines the IUPAC name of the product from the reaction
# of methyl phenyl sulfoxide with triflic anhydride and trimethylsilyl cyanide.
# The reaction is a Pummerer reaction.

# Step 1: Identify the product's structure.
# The Pummerer reaction of methyl phenyl sulfoxide (Ph-S(=O)-Me) results in
# an alpha-functionalized sulfide.
# The nucleophile is cyanide (CN-) from TMSCN.
# The product is Ph-S-CH2-CN.

# Step 2: Determine the IUPAC name.
# According to IUPAC rules, the nitrile group has higher priority than the sulfide.
# The parent compound is acetonitrile (CH3CN).
# The substituent at position 2 is the phenylsulfanyl group (Ph-S-).
# Ph = Phenyl
# S = Sulfa
# CH2-CN is derived from acetonitrile.

position = "2"
substituent = "Phenylsulfanyl"
parent = "acetonitrile"

# Constructing the IUPAC name
# The name format is position-(substituent)parent
final_iupac_name = f"{position}-({substituent}){parent}"

# Print the final IUPAC name. The numbers used in the name are also printed.
print(f"The number in the name is: {position}")
print(f"The final IUPAC name of the product is: {final_iupac_name}")
