import math

# 1. Define the input values for Chromium and Nickel equivalents.
cr_eq = 39
ni_eq = 29

# 2. Define the constants for the empirical formula.
# This formula is an approximation of the WRC-1992 diagram.
c1 = 3.34
c2 = 2.46
c3 = 28.6

# 3. Calculate the ferrite level using the formula.
ferrite_level = c1 * cr_eq - c2 * ni_eq - c3

# 4. Round the result to the nearest 10.
rounded_ferrite_level = int(round(ferrite_level / 10.0)) * 10

# 5. Print the calculation steps and the final result.
print("To find the approximate ferrite level, we use the formula: Ferrite Level = 3.34 * Cr_eq - 2.46 * Ni_eq - 28.6")
print("\nSubstituting the given values:")
print(f"{c1} * {cr_eq} - {c2} * {ni_eq} - {c3} = {ferrite_level:.2f}")
print(f"\nThe calculated approximate ferrite level is {ferrite_level:.2f}.")
print(f"Rounded to the nearest 10, the final answer is {rounded_ferrite_level}.")
