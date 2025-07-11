import math

# Part 1: How many pearls were there altogether?

# Number of pearls that remained on the string
# "a seven shy of eleven times eleven"
remained_on_string = (11 * 11) - 7

# The equation is: Total = (Total/6 + Total/5 + Total/3 + Total/10) + remained_on_string
# Rearranging to solve for Total:
# Total - (Total/6 + Total/5 + Total/3 + Total/10) = remained_on_string
# Total * (1 - (1/6 + 1/5 + 1/3 + 1/10)) = remained_on_string
# Total = remained_on_string / (1 - (1/6 + 1/5 + 1/3 + 1/10))

# Calculate the sum of the fractions
sum_of_fractions = 1/6 + 1/5 + 1/3 + 1/10

# Solve for the total number of pearls
# Using floating-point arithmetic and rounding to handle any precision issues
total_pearls = round(remained_on_string / (1 - sum_of_fractions))

print("--- Solving the Riddle ---")
print("First, let's find the total number of pearls in the necklace.")
print(f"The number of pearls remaining on the string is (11 * 11) - 7 = {remained_on_string}")

# Define the components of the scattered pearls
p_floor = total_pearls / 6
p_bed = total_pearls / 5
p_woman = total_pearls / 3
p_lover = total_pearls / 10

# Print the full equation with the solved values
print("\nThe final equation based on the riddle is:")
print(f"{p_floor:.0f} (fell to the floor) + {p_bed:.0f} (fell on the bed) + {p_woman:.0f} (saved by woman) + {p_lover:.0f} (caught by lover) + {remained_on_string} (remained on string) = {total_pearls}")

print(f"\nSo, there were {total_pearls} pearls altogether.")
print("\n--------------------------------")

# Part 2: How many more pearls are needed?

# Calculate the number of pearls that had fallen
fallen_pearls = total_pearls - remained_on_string

# They find 1/3rd of the fallen pearls
found_pearls = fallen_pearls / 3

# Calculate how many pearls they have now (on string + found)
current_pearl_count = remained_on_string + found_pearls

# Calculate how many more pearls are needed to get back to the original total
needed_pearls = total_pearls - current_pearl_count

print("\n--- The Second Question ---")
print(f"Total fallen pearls: {fallen_pearls:.0f}")
print(f"Pearls they find back (1/3 of fallen): {fallen_pearls:.0f} / 3 = {found_pearls:.0f}")
print(f"Total pearls they have now: {remained_on_string} (on string) + {found_pearls:.0f} (found) = {current_pearl_count:.0f}")
print(f"To restore the necklace, they need: {total_pearls} - {current_pearl_count:.0f} = {needed_pearls:.0f}")

print(f"\nThey will need {needed_pearls:.0f} more pearls.")

# The final answer for the second part of the question
final_answer = int(needed_pearls)
<<<304>>>