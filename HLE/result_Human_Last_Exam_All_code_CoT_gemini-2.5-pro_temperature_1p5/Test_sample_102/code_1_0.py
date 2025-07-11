import math

# Define the input values from the user prompt
chromium_equivalent = 39
nickel_equivalent = 29

# Define the coefficients for the WRC-1992 approximation formula
# Formula: FN = coeff_cr * Cr_eq - coeff_ni * Ni_eq - intercept
coeff_cr = 3.34
coeff_ni = 2.46
intercept = 28.6

# Calculate the ferrite number (FN)
ferrite_number = coeff_cr * chromium_equivalent - coeff_ni * nickel_equivalent - intercept

# Round the result to the nearest 10
# We can do this by dividing by 10, rounding to the nearest integer, and then multiplying by 10.
rounded_ferrite_level = int(round(ferrite_number / 10.0)) * 10

# Print the final equation with the rounded result, showing all numbers involved
print(f"The approximate ferrite level is calculated using the formula: FN = {coeff_cr} * Cr_eq - {coeff_ni} * Ni_eq - {intercept}")
print("Plugging in the values:")
print(f"{coeff_cr} * {chromium_equivalent} - {coeff_ni} * {nickel_equivalent} - {intercept} = {ferrite_number:.2f}")
print(f"Rounding to the nearest 10 gives: {rounded_ferrite_level}")
print("\nFinal Equation and Answer:")
print(f"{coeff_cr} * {chromium_equivalent} - {coeff_ni} * {nickel_equivalent} - {intercept} = {rounded_ferrite_level}")
<<<30>>>