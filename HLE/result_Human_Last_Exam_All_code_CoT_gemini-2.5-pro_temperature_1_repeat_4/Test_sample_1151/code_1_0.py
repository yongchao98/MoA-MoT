import math

# Step 1: Calculate the number of pearls remaining on the string.
# "a seven shy of eleven times eleven"
remained_on_string = (11 * 11) - 7

# Step 2: Set up the equation to find the total number of pearls (x).
# x = (1/6)x + (1/5)x + (1/3)x + (1/10)x + remained_on_string
# x - (1/6 + 1/5 + 1/3 + 1/10)x = remained_on_string
# x * (1 - (1/6 + 1/5 + 1/3 + 1/10)) = remained_on_string
# x = remained_on_string / (1 - (1/6 + 1/5 + 1/3 + 1/10))

# The sum of the fractions of pearls that fell
fraction_fallen = (1/6) + (1/5) + (1/3) + (1/10)

# Solve for x, the total number of pearls
# Using floating point numbers and rounding to the nearest integer
total_pearls = remained_on_string / (1 - fraction_fallen)
total_pearls = int(round(total_pearls))

print(f"How many pearls were there altogether?")
print(f"Answer: {total_pearls}\n")

# Step 3: Display the equation with the final numbers.
floor_pearls = int(round(total_pearls / 6))
bed_pearls = int(round(total_pearls / 5))
woman_saved = int(round(total_pearls / 3))
lover_caught = int(round(total_pearls / 10))

print("The final equation with all numbers is:")
print(f"{total_pearls} = {floor_pearls} (floor) + {bed_pearls} (bed) + {woman_saved} (woman) + {lover_caught} (lover) + {remained_on_string} (string)\n")

# Step 4: Calculate the number of pearls needed for the second question.
total_fallen = total_pearls - remained_on_string
found_back = math.floor(total_fallen / 3)
still_needed = total_fallen - found_back

print("How many more are they now gonna need if they manage to find back only 1/3rd of the fallen ones?")
print(f"Answer: {still_needed}")
