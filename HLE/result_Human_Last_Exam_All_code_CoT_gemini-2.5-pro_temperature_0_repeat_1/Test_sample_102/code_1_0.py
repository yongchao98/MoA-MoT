import math

# Define the nickel and chromium equivalents
ni_eq = 29
cr_eq = 39

# Calculate the ferrite level using the AWS formula
# Formula: Ferrite Level = 3.34 * Cr_eq - 2.46 * Ni_eq - 28.6
ferrite_level = 3.34 * cr_eq - 2.46 * ni_eq - 28.6

# Round the result to the nearest 10
# To round to the nearest 10, we divide by 10, round to the nearest integer, and then multiply by 10.
rounded_ferrite_level = round(ferrite_level / 10) * 10

# Print the calculation steps and the final result
print(f"To find the ferrite level, we use the formula: Ferrite Level = 3.34 * Cr_eq - 2.46 * Ni_eq - 28.6")
print(f"Plugging in the values Cr_eq = {cr_eq} and Ni_eq = {ni_eq}:")
print(f"Ferrite Level = 3.34 * {cr_eq} - 2.46 * {ni_eq} - 28.6")
print(f"Calculated Ferrite Level = {ferrite_level:.2f}")
print(f"The approximate ferrite level rounded to the nearest 10 is: {rounded_ferrite_level}")