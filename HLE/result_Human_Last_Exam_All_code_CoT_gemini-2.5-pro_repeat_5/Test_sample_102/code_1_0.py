import math

# Step 1: Define the input values
ni_eq = 29.0
cr_eq = 39.0

# Step 2: Use an empirical formula to calculate the ferrite level (FN)
# FN = 4.47 * Cr_eq - 3.44 * Ni_eq - 28.52
term1 = 4.47 * cr_eq
term2 = 3.44 * ni_eq
constant = 28.52

# Step 3: Calculate the result
ferrite_level = term1 - term2 - constant

# Step 4: Print the full equation and the result
print(f"Using the formula: FN = 4.47 * Cr_eq - 3.44 * Ni_eq - 28.52")
print(f"Calculation: {4.47} * {cr_eq} - {3.44} * {ni_eq} - {constant} = {ferrite_level:.2f}")

# Step 5: Round the result to the nearest 10
rounded_ferrite = round(ferrite_level / 10) * 10

print(f"\nThe approximate ferrite level is {ferrite_level:.2f}%.")
print(f"Rounded to the nearest 10, the ferrite level is {rounded_ferrite}%.")
