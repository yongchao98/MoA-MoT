import math

# --- Part 1: How many pearls were there altogether? ---

# Calculate the number of pearls that remained on the string
# "a seven shy of eleven times eleven"
eleven_squared = 11 * 11
remaining_on_string = eleven_squared - 7

# The fractions of the total pearls (x) that fell
# 1/6 + 1/5 + 1/3 + 1/10
fraction_fallen = 1/6 + 1/5 + 1/3 + 1/10

# The fraction that remained on the string is 1 minus the fraction that fell
fraction_remaining = 1 - fraction_fallen

# The total number of pearls is the number remaining divided by the fraction remaining
# fraction_remaining * total_pearls = remaining_on_string
total_pearls = remaining_on_string / fraction_remaining
# Ensure the result is an integer, as we are counting pearls
total_pearls = int(round(total_pearls))

print("--- Part 1: How many pearls were there altogether? ---")
print("Let 'x' be the total number of pearls.")
print("The equation is based on the parts of the necklace:")
print("x = (pearls on floor) + (pearls on bed) + (saved by woman) + (caught by lover) + (remaining on string)")
# The problem asks to output each number in the final equation.
print(f"x = (1/6)x + (1/5)x + (1/3)x + (1/10)x + ({eleven_squared} - 7)")
print(f"x = ({fraction_fallen:.4f})x + {remaining_on_string}")
print(f"(1 - {fraction_fallen:.4f})x = {remaining_on_string}")
print(f"({fraction_remaining:.4f})x = {remaining_on_string}")
print(f"x = {remaining_on_string} / {fraction_remaining:.4f}")
print(f"\nThere were {total_pearls} pearls altogether.\n")


# --- Part 2: How many more pearls are needed? ---

# Calculate the number of pearls that had fallen
fallen_pearls = total_pearls - remaining_on_string

# Calculate the number of pearls they found back (1/3rd of fallen)
found_pearls = fallen_pearls / 3
# Use math.floor since you can't find a fraction of a pearl
found_pearls = math.floor(found_pearls)

# Calculate how many pearls they have in total now
current_pearls_count = remaining_on_string + found_pearls

# Calculate how many more they need to complete the necklace
needed_pearls = total_pearls - current_pearls_count

print("--- Part 2: How many more pearls are needed? ---")
print(f"Number of fallen pearls: {total_pearls} - {remaining_on_string} = {fallen_pearls}")
print(f"Number of pearls they found back (1/3 of fallen): {fallen_pearls} / 3 = {found_pearls}")
print(f"Total pearls they have now: {remaining_on_string} (on string) + {found_pearls} (found) = {current_pearls_count}")
print(f"Number of pearls still needed: {total_pearls} (original) - {current_pearls_count} (current) = {needed_pearls}")
print(f"\nThey will need {needed_pearls} more pearls.")

<<<304>>>