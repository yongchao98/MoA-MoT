import math

# Step 1: Calculate the number of pearls remaining on the string.
# "a seven shy of eleven times eleven"
remaining_on_string = (11 * 11) - 7

# Step 2 & 3: Solve for the total number of pearls (x).
# The equation is: x = (x/6) + (x/5) + (x/3) + (x/10) + remaining_on_string
# This can be rewritten to isolate x:
# x - (x/6 + x/5 + x/3 + x/10) = remaining_on_string
# To simplify, we find a common denominator for the fractions (30):
# x - (5x/30 + 6x/30 + 10x/30 + 3x/30) = remaining_on_string
# x - (24x/30) = remaining_on_string
# (6x/30) = remaining_on_string -> (x/5) = remaining_on_string
# x = remaining_on_string * 5

total_pearls = remaining_on_string * 5

# Calculate each part for the final equation printout
floor_pearls = total_pearls / 6
bed_pearls = total_pearls / 5
woman_saved_pearls = total_pearls / 3
lover_caught_pearls = total_pearls / 10

# Print the solution to the first question
print("To find the total number of pearls, we solve the following equation based on the riddle:")
# Using math.ceil to handle any potential floating point inaccuracies and ensure integer output
print(f"{math.ceil(total_pearls)} = {math.ceil(floor_pearls)} (floor) + {math.ceil(bed_pearls)} (bed) + {math.ceil(woman_saved_pearls)} (woman) + {math.ceil(lover_caught_pearls)} (lover) + {math.ceil(remaining_on_string)} (string)")
print(f"\nThere were {math.ceil(total_pearls)} pearls on the necklace altogether.\n")


# Step 4: Calculate the total number of fallen pearls
fallen_pearls = total_pearls - remaining_on_string

# Step 5: Calculate how many of the fallen pearls were found
found_pearls = fallen_pearls / 3

# Step 6: Calculate how many pearls are still needed
needed_pearls = fallen_pearls - found_pearls

# Print the solution to the second question
print("For the second part of the question:")
print(f"Number of fallen pearls: {math.ceil(fallen_pearls)}")
print(f"Number of pearls they find back (1/3rd): {math.ceil(found_pearls)}")
print(f"Thus, the number of pearls they still need is {math.ceil(fallen_pearls)} - {math.ceil(found_pearls)} = {math.ceil(needed_pearls)}.")

<<<304>>>