# --- FMO Theory Prediction for Competing Electrocyclizations ---
# This script predicts the product ratio for the thermal reaction of
# (2Z,4Z,6Z,8E)-deca-2,4,6,8-tetraene.

# The reaction has two competing pathways:
# 1. An 8-pi electron conrotatory cyclization, leading to cis-isomer A.
# 2. A 6-pi electron disrotatory cyclization, leading to trans-isomer B.

# Number of pi electrons involved in the formation of product A
pi_electrons_A = 8

# Number of pi electrons involved in the formation of product B
pi_electrons_B = 6

# The ratio of products A/B is determined by the relative rates of the
# two competing reactions. Generally, reactions involving fewer pi electrons
# are kinetically favored. We will model the ratio of the products
# as the inverse ratio of the pi electrons involved.
# Ratio A/B = k_A / k_B ≈ pi_electrons_B / pi_electrons_A

# Calculate the ratio
ratio_A_to_B = pi_electrons_B / pi_electrons_A

# Print the reasoning and the final calculation
print("Prediction of the product ratio A/B using Frontier Molecular Orbital theory.")
print("Product A is formed via an 8-pi electron electrocyclization.")
print("Product B is formed via a 6-pi electron electrocyclization.")
print("The kinetic ratio of A/B can be estimated as the ratio of pi electrons involved in the competing pathways.")
print(f"Ratio A/B = (pi electrons for B) / (pi electrons for A)")
print(f"Ratio A/B = {pi_electrons_B} / {pi_electrons_A}")
print(f"The calculated ratio of A to B is: {ratio_A_to_B}")
