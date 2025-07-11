import math

# Step 1: Define the input values for chromium and nickel equivalents.
cr_eq = 39.0
ni_eq = 29.0

# Step 2: Apply Suutala's empirical formula to calculate the ferrite level.
# The formula is: Ferrite Level = 3.34 * Cr_eq - 2.86 * Ni_eq - 27.5
ferrite_level = 3.34 * cr_eq - 2.86 * ni_eq - 27.5

# Step 3: Round the calculated ferrite level to the nearest 10.
# This is done by dividing by 10, rounding to the nearest whole number, and then multiplying by 10.
rounded_ferrite_level = int(round(ferrite_level / 10.0)) * 10

# Step 4: Print the full equation and the final answer.
# The following lines will display each number used in the equation and the final rounded result.
print(f"Calculating using the formula: 3.34 * Cr_eq - 2.86 * Ni_eq - 27.5")
print(f"The equation with the given values is: 3.34 * {cr_eq} - 2.86 * {ni_eq} - 27.5 = {ferrite_level:.2f}")
print(f"The approximate ferrite level is {ferrite_level:.2f}%.")
print(f"Rounding to the nearest 10, the final ferrite level is: {rounded_ferrite_level}")