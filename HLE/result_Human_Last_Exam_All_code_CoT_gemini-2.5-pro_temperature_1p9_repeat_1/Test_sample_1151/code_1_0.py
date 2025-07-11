import math

# Part 1: How many pearls were there altogether?

# Calculate the number of pearls remaining on the string
pearls_remaining_on_string = 11 * 11 - 7

# The sum of the fractions of the total pearls that fell
sum_of_fractions = 1/6 + 1/5 + 1/3 + 1/10

# The equation is: total_pearls = (sum_of_fractions * total_pearls) + pearls_remaining_on_string
# We can rearrange this to solve for total_pearls:
# total_pearls * (1 - sum_of_fractions) = pearls_remaining_on_string
# total_pearls = pearls_remaining_on_string / (1 - sum_of_fractions)
total_pearls = pearls_remaining_on_string / (1 - sum_of_fractions)

# Since the number of pearls must be a whole number, we cast it to an integer.
total_pearls = int(round(total_pearls))

print("--- Part 1: Total Pearls on the Necklace ---")
# Calculate the number of pearls for each fraction
p_floor = total_pearls / 6
p_bed = total_pearls / 5
p_woman = total_pearls / 3
p_lover = total_pearls / 10

# Print the final equation with all numbers filled in
print("The equation representing the story is:")
print(f"{total_pearls} = {int(p_floor)} (on floor) + {int(p_bed)} (on bed) + {int(p_woman)} (saved by woman) + {int(p_lover)} (caught by lover) + {pearls_remaining_on_string} (on string)")

print(f"\nTherefore, there were {total_pearls} pearls altogether on the necklace.")
print("-" * 50)


# Part 2: How many more pearls are needed?
target_pearls = 500

# Calculate the total number of fallen pearls
fallen_pearls = total_pearls - pearls_remaining_on_string

# They find 1/3rd of the fallen ones
found_pearls = fallen_pearls / 3

# The total number of pearls they have now
current_pearls = pearls_remaining_on_string + found_pearls

# The number of additional pearls needed to reach their target of 500
needed_pearls = target_pearls - current_pearls

print("--- Part 2: Pearls Needed for a New Necklace ---")
print(f"The lovers had {int(fallen_pearls)} fallen pearls and found 1/3rd of them, which is {int(found_pearls)} pearls.")
print(f"They now have {int(pearls_remaining_on_string)} (on the string) + {int(found_pearls)} (found) = {int(current_pearls)} pearls.")
print(f"To make their desired {target_pearls}-pearl necklace, they will need {int(needed_pearls)} more pearls.")
