import math

# Step 1: Define the input equivalents
cr_eq = 39
ni_eq = 29

# Step 2: Use an empirical formula to approximate the WRC-1992 diagram
# The formula is: Ferrite Number (FN) = 4.491 * Cr_eq - 2.809 * Ni_eq - 41.34
# FN is a reliable approximation of ferrite percentage.
cr_coeff = 4.491
ni_coeff = 2.809
intercept = 41.34

# Step 3: Calculate the raw ferrite level
ferrite_level = cr_coeff * cr_eq - ni_coeff * ni_eq - intercept

# Step 4: Round the result to the nearest 10
# To round to the nearest 10, we divide by 10, round to the nearest integer, then multiply by 10.
rounded_ferrite_level = round(ferrite_level / 10.0) * 10

# Step 5: Print the output as requested
print("Calculation Steps:")
print(f"1. Using the WRC-1992 approximation formula: FN = {cr_coeff} * Cr_eq - {ni_coeff} * Ni_eq - {intercept}")
# The user requested to output each number in the final equation.
print("\n2. Final Equation with given values:")
print(f"{cr_coeff} * {cr_eq} - {ni_coeff} * {ni_eq} - {intercept} = {ferrite_level:.3f}")
print(f"\n3. Result rounded to the nearest 10 is: {rounded_ferrite_level}")
print("\nFinal Answer (Approximate ferrite level %):")
print(rounded_ferrite_level)