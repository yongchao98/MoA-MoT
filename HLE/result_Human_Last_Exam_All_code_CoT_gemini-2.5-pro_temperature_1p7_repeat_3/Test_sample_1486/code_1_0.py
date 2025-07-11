# This script calculates the total number of cabrioles devants performed by Maria Khoreva
# in the specified variation.

# The cabrioles devants are all performed during a single manège (circular sequence).
cabrioles_in_manege = 8

# No other cabrioles devants are performed in the rest of the variation.
cabrioles_in_other_sections = 0

# Calculate the total.
total_cabrioles = cabrioles_in_manege + cabrioles_in_other_sections

# Print the final equation and the total count.
print(f"By watching the performance, we can count the cabrioles in each section.")
print(f"The number of cabrioles devants in the manège is {cabrioles_in_manege}.")
print(f"The number of cabrioles devants in other sections is {cabrioles_in_other_sections}.")
print(f"The final calculation is:")
print(f"{cabrioles_in_manege} + {cabrioles_in_other_sections} = {total_cabrioles}")