import math

# Step 1: Define the input values
ni_eq = 29.0
cr_eq = 39.0

# Step 2: Apply the WRC-1992 approximation formula
# The formula is: Ferrite Level = 3.34 * Cr_eq - 2.46 * Ni_eq - 28.6
term1 = 3.34 * cr_eq
term2 = 2.46 * ni_eq
constant = 28.6
ferrite_level = term1 - term2 - constant

# Step 3: Round the result to the nearest 10
# We can do this by dividing by 10, rounding to the nearest integer, and then multiplying by 10.
rounded_ferrite_level = int(round(ferrite_level / 10.0)) * 10

# Step 4: Print the equation and the final answer
print("Calculating the approximate ferrite level using the WRC-1992 diagram approximation.")
print("Formula: Ferrite Level = 3.34 * Cr_eq - 2.46 * Ni_eq - 28.6")
print(f"Substituting the given values:")
print(f"Ferrite Level = 3.34 * {cr_eq} - 2.46 * {ni_eq} - {constant}")
print(f"Calculated Value = {ferrite_level:.2f}")
print(f"Rounding to the nearest 10, the approximate ferrite level is: {rounded_ferrite_level}")
<<<30>>>