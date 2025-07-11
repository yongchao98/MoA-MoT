import math

# Step 1 & 2: Set up the equation from the riddle's details.
# Let 'x' be the total number of pearls.

# The number of pearls that remained on the string is "seven shy of eleven times eleven".
eleven_squared = 11 * 11
remained_on_string = eleven_squared - 7

# The fractions of the total pearls that fell.
# 1/6 to the floor, 1/5 to the bed, 1/3 saved by the woman, 1/10 caught by the lover.
# The equation is: x = (1/6)x + (1/5)x + (1/3)x + (1/10)x + remained_on_string
# We can find x by first isolating it:
# x - (1/6 + 1/5 + 1/3 + 1/10)x = remained_on_string
# To sum the fractions, we find a common denominator, which is 30.
# (5/30 + 6/30 + 10/30 + 3/30)x = (24/30)x
# So, x - (24/30)x = remained_on_string
# (6/30)x = remained_on_string, which simplifies to (1/5)x = remained_on_string
# Therefore, x = remained_on_string * 5

# Step 3: Solve for 'x', the total number of pearls.
total_pearls = remained_on_string * 5

print("--- Part 1: How many pearls were there altogether? ---")
print("The riddle can be expressed as a single equation. Let 'x' be the total number of pearls.")
print("The final equation is: x = (x / 6) + (x / 5) + (x / 3) + (x / 10) + (11*11 - 7)")
print(f"The number of pearls remaining on the string is 11*11 - 7 = {remained_on_string}.")
print(f"Solving the equation for x gives the total number of pearls.")
print(f"The total number of pearls on the necklace was: {total_pearls}\n")


# Step 4: Answer the second question.
# How many more do they need if they find 1/3rd of the fallen ones?

# First, calculate the number of fallen pearls.
fallen_pearls = total_pearls - remained_on_string

# Next, calculate how many they found.
found_pearls = fallen_pearls / 3

# Finally, calculate how many more they need, which is the number of pearls still lost.
needed_pearls = fallen_pearls - found_pearls

print("--- Part 2: How many more pearls do they need? ---")
print(f"First, calculate the number of pearls that fell off: {total_pearls} (total) - {remained_on_string} (on string) = {int(fallen_pearls)} pearls.")
print(f"They manage to find 1/3rd of the fallen pearls: {int(fallen_pearls)} / 3 = {int(found_pearls)} pearls.")
print(f"To restore the necklace, they still need to find the rest of the fallen pearls.")
print(f"The number of pearls they still need is: {int(fallen_pearls)} - {int(found_pearls)} = {int(needed_pearls)}.")

# The final answer in the required format.
# It should contain the total number of pearls and the number of pearls they still need.
final_answer = f"<<<{total_pearls}, {int(needed_pearls)}>>>"