import math

# --- Part 1: How many pearls were there altogether? ---

# First, calculate the number of pearls that remained on the string.
# The riddle states "a seven shy of eleven times eleven".
eleven = 11
seven = 7
remained_on_string = (eleven * eleven) - seven

# The total number of pearls, X, can be found with the equation:
# X = (1/6)X + (1/5)X + (1/3)X + (1/10)X + remained_on_string
#
# To solve for X, we rearrange the formula:
# X - (1/6)X - (1/5)X - (1/3)X - (1/10)X = remained_on_string
# X * (1 - 1/6 - 1/5 - 1/3 - 1/10) = remained_on_string
# X = remained_on_string / (1 - 1/6 - 1/5 - 1/3 - 1/10)

# Calculate the fractional coefficient of X
coefficient = 1 - (1/6 + 1/5 + 1/3 + 1/10)

# Calculate the total number of pearls
# We round to the nearest integer as the number of pearls must be a whole number.
total_pearls = int(round(remained_on_string / coefficient))

# Calculate the number of pearls for each fraction mentioned in the riddle
p_floor = total_pearls / 6
p_bed = total_pearls / 5
p_woman = total_pearls / 3
p_lover = total_pearls / 10

print("--- How many pearls were there altogether? ---")
print(f"The number of pearls that remained on the string was {eleven} * {eleven} - {seven} = {remained_on_string}")
print(f"The total number of pearls on the necklace was: {total_pearls}")

print("\nThis is based on the final equation:")
# We cast the numbers to int for clean printing
print(f"{int(p_floor)} (fell to floor) + {int(p_bed)} (fell on bed) + {int(p_woman)} (saved by woman) + {int(p_lover)} (caught by lover) + {remained_on_string} (remained on string) = {total_pearls}")

# --- Part 2: How many more pearls are needed? ---

# Calculate the total number of pearls that had fallen
total_fallen = total_pearls - remained_on_string

# Calculate how many of the fallen pearls they find back (1/3rd)
found_pearls = total_fallen / 3

# Calculate the total number of pearls they have now
current_pearls_total = remained_on_string + found_pearls

# Calculate how many more pearls they need to get back to the original total
needed_pearls = total_pearls - current_pearls_total

print("\n--- How many more pearls are they now gonna need? ---")
print(f"Of the {total_pearls} pearls, {int(total_fallen)} had fallen.")
print(f"They manage to find back 1/3 of the fallen pearls: {int(total_fallen)} / 3 = {int(found_pearls)}")
print(f"The number of additional pearls needed is: {int(needed_pearls)}")
