import math

# Step 1: Define the given Nickel and Chromium equivalent values.
ni_eq = 29
cr_eq = 39

# Step 2: Define the coefficients for the Schaeffler diagram approximation formula.
# The formula is: Ferrite % = a * Cr_eq - b * Ni_eq - c
a = 4.1
b = 2.9
c = 35

# Step 3: Print the full equation with the values substituted.
print(f"Using the formula: Ferrite % = {a} * Cr_eq - {b} * Ni_eq - {c}")
print(f"Plugging in the values: Ferrite % = {a} * {cr_eq} - {b} * {ni_eq} - {c}")

# Step 4: Calculate the ferrite percentage.
ferrite_percent = a * cr_eq - b * ni_eq - c

# Step 5: Round the result to the nearest 10.
# We do this by dividing by 10, rounding to the nearest integer, and then multiplying by 10.
rounded_ferrite = round(ferrite_percent / 10) * 10

# Step 6: Print the final answer.
print(f"The calculated ferrite level is approximately {ferrite_percent:.1f}%.")
print(f"Rounded to the nearest 10, the ferrite level is {rounded_ferrite}%.")
