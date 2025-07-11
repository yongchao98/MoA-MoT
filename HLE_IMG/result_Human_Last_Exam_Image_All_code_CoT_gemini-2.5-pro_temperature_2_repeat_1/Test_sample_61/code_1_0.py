# This script calculates the molecular formula of product A (2-benzyl-1-indanone).

# The structure of 2-benzyl-1-indanone can be broken down into:
# 1. An indanone core skeleton.
# 2. A benzyl group attached at position 2.
# Let's count the atoms for each part and sum them up.

# 1. Count Carbon atoms (C)
# The indanone core consists of a C6H4 benzene ring fused with a 5-membered ring.
# The five-membered ring has a carbonyl group (C=O), an alpha-carbon (CH), and a beta-carbon (CH2).
c_benzene_ring = 6
c_carbonyl = 1
c_alpha_carbon = 1
c_beta_carbon = 1
# The benzyl group (C7H7) has 7 carbon atoms.
c_benzyl_group = 7
# Total carbon atoms
total_c = c_benzene_ring + c_carbonyl + c_alpha_carbon + c_beta_carbon + c_benzyl_group
print("Calculation for Carbon atoms:")
print(f"{c_benzene_ring} (from indanone's benzene ring) + {c_carbonyl} (from C=O group) + {c_alpha_carbon} (from alpha-carbon) + {c_beta_carbon} (from beta-carbon) + {c_benzyl_group} (from benzyl group) = {total_c}")

# 2. Count Hydrogen atoms (H)
# In the C6H4 part of the indanone ring:
h_benzene_ring = 4
# On the alpha-carbon (CH):
h_alpha_carbon = 1
# On the beta-carbon (CH2):
h_beta_carbon = 2
# In the benzyl group (C7H7):
h_benzyl_group = 7
# Total hydrogen atoms
total_h = h_benzene_ring + h_alpha_carbon + h_beta_carbon + h_benzyl_group
print("\nCalculation for Hydrogen atoms:")
print(f"{h_benzene_ring} (from C6H4 ring) + {h_alpha_carbon} (from alpha-carbon) + {h_beta_carbon} (from beta-carbon) + {h_benzyl_group} (from benzyl group) = {total_h}")

# 3. Count Oxygen atoms (O)
# There is one ketone group (C=O) in the indanone core.
total_o = 1
print("\nCalculation for Oxygen atoms:")
print(f"1 (from C=O group) = {total_o}")

# 4. Final Molecular Formula
print(f"\nBased on the calculations, the molecular formula of compound A is C{total_c}H{total_h}O{total_o}.")