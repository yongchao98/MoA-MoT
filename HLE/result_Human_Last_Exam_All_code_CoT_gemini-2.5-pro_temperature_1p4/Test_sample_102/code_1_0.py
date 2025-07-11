import math

# Step 1: Define the given equivalent values.
cr_eq = 39
ni_eq = 29

# Step 2 & 3: Use the Schaeffler diagram-based formula to calculate the ferrite percentage.
# The formula is: Ferrite % = 4 * Cr_eq - 2.8 * Ni_eq - 30
term1_val = 4 * cr_eq
term2_val = 2.8 * ni_eq
constant = 30

# Step 4: Calculate the result.
ferrite_percent = term1_val - term2_val - constant

# Step 5: Round the result to the nearest 10.
# We can achieve this by dividing by 10, rounding to the nearest integer, and then multiplying by 10.
rounded_ferrite_level = int(round(ferrite_percent / 10.0)) * 10

# Step 6: Print the output in the required format.
print("The calculation for the ferrite level using the formula is:")
print(f"4 * {cr_eq} - 2.8 * {ni_eq} - {constant} = {ferrite_percent:.1f}")
print("\nThe approximate ferrite level, rounded to the nearest 10, is:")
print(rounded_ferrite_level)