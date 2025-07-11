import math

# Step 1: Define the given equivalent values
cr_eq = 39
ni_eq = 29

# Step 2: Use an empirical formula to approximate the ferrite percentage.
# This formula is a linear approximation of the WRC-1992 diagram.
# Formula: Ferrite % = 3.34 * Cr_eq - 2.46 * Ni_eq - 28.6
a = 3.34
b = 2.46
c = 28.6

# Step 3: Calculate the ferrite percentage
ferrite_percentage = a * cr_eq - b * ni_eq - c

# Step 4: Round the result to the nearest 10
# This is done by dividing by 10, rounding to the nearest whole number, and multiplying by 10.
rounded_ferrite_level = int(round(ferrite_percentage / 10.0)) * 10

# Step 5: Print the full equation and the final answer
print("Calculating the approximate ferrite level using the formula:")
print(f"Ferrite % = {a} * {cr_eq} - {b} * {ni_eq} - {c}")

print(f"\nCalculated result before rounding: {ferrite_percentage:.2f}%")
print(f"Final approximate ferrite level rounded to the nearest 10 is: {rounded_ferrite_level}")