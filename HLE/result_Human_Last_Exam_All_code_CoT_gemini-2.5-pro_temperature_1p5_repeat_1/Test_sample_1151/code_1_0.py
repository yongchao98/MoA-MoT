import math

# Part 1: How many pearls were there altogether?

# First, calculate the number of pearls remaining on the string
# "a seven shy of eleven times eleven"
remaining_on_string = (11 * 11) - 7

# The total number of pearls 'x' is the sum of its parts.
# x = (1/6)x + (1/5)x + (1/3)x + (1/10)x + remaining_on_string
# To solve for x, we can rearrange the equation:
# x - (1/6)x - (1/5)x - (1/3)x - (1/10)x = remaining_on_string
# x * (1 - 1/6 - 1/5 - 1/3 - 1/10) = remaining_on_string
# x = remaining_on_string / (1 - (1/6 + 1/5 + 1/3 + 1/10))

# The sum of the fractions of fallen pearls
sum_of_fractions = (1/6) + (1/5) + (1/3) + (1/10)

# Solve for the total number of pearls
# Using floating-point numbers to ensure precision, then converting to int
total_pearls = remaining_on_string / (1 - sum_of_fractions)
total_pearls = int(round(total_pearls))

# Calculate the number of pearls for each part based on the total
floor_pearls = int(total_pearls * (1/6))
bed_pearls = int(total_pearls * (1/5))
woman_saved_pearls = int(total_pearls * (1/3))
lover_caught_pearls = int(total_pearls * (1/10))

print("--- The Riddle of the Pearls ---")
print("\nSolving for the total number of pearls (x):")
print(f"The equation is: x = (1/6)x + (1/5)x + (1/3)x + (1/10)x + ((11*11)-7)")
print(f"\nBreaking down the necklace:")
print(f"Pearls that fell to the floor: {floor_pearls}")
print(f"Pearls that fell on the bed: {bed_pearls}")
print(f"Pearls the woman saved: {woman_saved_pearls}")
print(f"Pearls the lover caught: {lover_caught_pearls}")
print(f"Pearls that remained on the string: {remaining_on_string}")

print("\nThe final equation with the numbers is:")
print(f"{floor_pearls} + {bed_pearls} + {woman_saved_pearls} + {lover_caught_pearls} + {remaining_on_string} = {total_pearls}")

print(f"\nAnswer 1: There were {total_pearls} pearls altogether.")

# Part 2: How many more pearls do they need?
fallen_pearls = total_pearls - remaining_on_string
found_pearls = math.floor(fallen_pearls / 3)
pearls_still_needed = total_pearls - (remaining_on_string + found_pearls)

print("\n--- The Second Question ---")
print(f"Total pearls fallen from the string: {fallen_pearls}")
print(f"Number of fallen pearls they find back (1/3rd): {found_pearls}")
print(f"Total pearls they have now (remaining on string + found): {remaining_on_string + found_pearls}")

print(f"\nAnswer 2: They will need {pearls_still_needed} more pearls.")
