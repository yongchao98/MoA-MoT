import math

# Part 1: Finding the total number of pearls on the necklace.

# Calculate the number of pearls that remained on the string.
# "a seven shy of eleven times eleven"
pearls_on_string = (11 * 11) - 7

# Let 'x' be the total number of pearls. The equation is:
# x = (1/6)x + (1/5)x + (1/3)x + (1/10)x + pearls_on_string
# To solve for x, we rearrange the equation:
# x - (1/6 + 1/5 + 1/3 + 1/10)x = pearls_on_string
# x * (1 - (1/6 + 1/5 + 1/3 + 1/10)) = pearls_on_string
# x = pearls_on_string / (1 - (1/6 + 1/5 + 1/3 + 1/10))

# The sum of the fractions of the fallen pearls
sum_of_fractions = 1/6 + 1/5 + 1/3 + 1/10

# Calculate the total number of pearls
# We use round() to handle potential floating-point inaccuracies and get a whole number.
total_pearls = round(pearls_on_string / (1 - sum_of_fractions))

print("--- Part 1: How many pearls were there altogether? ---")

# Calculate the number for each part of the riddle
floor_pearls = math.ceil(total_pearls * (1/6))
bed_pearls = math.ceil(total_pearls * (1/5))
woman_pearls = math.ceil(total_pearls * (1/3))
lover_pearls = math.ceil(total_pearls * (1/10))

print(f"The final equation showing the distribution of pearls is:")
# Note: Using integer values for display clarity
print(f"{int(floor_pearls)} (on floor) + {int(bed_pearls)} (on bed) + {int(woman_pearls)} (woman saved) + {int(lover_pearls)} (lover caught) + {int(pearls_on_string)} (on string) = {int(total_pearls)}")
print(f"\nAnswer: There were {int(total_pearls)} pearls altogether.")

print("\n--- Part 2: How many more pearls are needed? ---")

# Calculate the number of pearls that fell off the string.
fallen_pearls = total_pearls - pearls_on_string
print(f"Total pearls that fell: {int(total_pearls)} - {int(pearls_on_string)} = {int(fallen_pearls)}")

# Calculate the number of pearls they found (1/3rd of the fallen ones).
found_pearls = fallen_pearls / 3
print(f"Pearls they found back: {int(fallen_pearls)} / 3 = {int(found_pearls)}")

# Calculate the total number of pearls they currently possess.
current_pearls_total = pearls_on_string + found_pearls
print(f"Current total pearls they have: {int(pearls_on_string)} + {int(found_pearls)} = {int(current_pearls_total)}")

# Calculate the number of pearls they still need.
needed_pearls = total_pearls - current_pearls_total
print(f"Pearls still needed to complete the necklace: {int(total_pearls)} - {int(current_pearls_total)} = {int(needed_pearls)}")
print(f"\nAnswer: They are going to need {int(needed_pearls)} more pearls.")

<<<304>>>