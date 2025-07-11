# Part 1: How many pearls were there altogether?

# The number of pearls remaining on the string
# "a seven shy of eleven times eleven"
remaining_on_string = (11 * 11) - 7

# The sum of the fractions of the scattered pearls
# 1/6 + 1/5 + 1/3 + 1/10
# The common denominator for 6, 5, 3, 10 is 30.
# 5/30 + 6/30 + 10/30 + 3/30 = 24/30 = 4/5
sum_of_fractions = 1/6 + 1/5 + 1/3 + 1/10

# The equation is: x = (sum_of_fractions * x) + remaining_on_string
# x = (4/5)x + 114
# x - (4/5)x = 114
# (1/5)x = 114
# x = 114 * 5
total_pearls = remaining_on_string / (1 - sum_of_fractions)
total_pearls = int(total_pearls)

print(f"First, we calculate the number of pearls remaining on the string: (11 * 11) - 7 = {remaining_on_string}")
print(f"The total number of pearls (x) can be found with the equation: x = (1/6)x + (1/5)x + (1/3)x + (1/10)x + {remaining_on_string}")
print(f"Solving for x, we find that the total number of pearls was {total_pearls}.")
print("\n--- Breakdown of the Pearls ---")

# Calculate the number of pearls for each fraction
floor_pearls = int(total_pearls / 6)
bed_pearls = int(total_pearls / 5)
woman_pearls = int(total_pearls / 3)
lover_pearls = int(total_pearls / 10)

print(f"On the floor (1/6): {floor_pearls}")
print(f"On the bed (1/5): {bed_pearls}")
print(f"Saved by the woman (1/3): {woman_pearls}")
print(f"Caught by the lover (1/10): {lover_pearls}")
print(f"Remaining on the string: {remaining_on_string}")

# Verify the total
print("\nThe final equation proving the total:")
print(f"{floor_pearls} + {bed_pearls} + {woman_pearls} + {lover_pearls} + {remaining_on_string} = {total_pearls}")

# Part 2: How many more are they now gonna need?

# Calculate the number of fallen pearls
fallen_pearls = total_pearls - remaining_on_string

# They find 1/3rd of the fallen ones
found_pearls = int(fallen_pearls / 3)

# Calculate how many pearls they have in total now
current_pearls = remaining_on_string + found_pearls

# Calculate how many more they need to complete the necklace
needed_pearls = total_pearls - current_pearls

print("\n--- Pearls Needed ---")
print(f"Total fallen pearls: {total_pearls} - {remaining_on_string} = {fallen_pearls}")
print(f"They find 1/3 of the fallen pearls: {fallen_pearls} / 3 = {found_pearls}")
print(f"Total pearls they have now (remaining + found): {remaining_on_string} + {found_pearls} = {current_pearls}")
print(f"To restore the necklace, they will need: {total_pearls} - {current_pearls} = {needed_pearls} more pearls.")

<<<304>>>