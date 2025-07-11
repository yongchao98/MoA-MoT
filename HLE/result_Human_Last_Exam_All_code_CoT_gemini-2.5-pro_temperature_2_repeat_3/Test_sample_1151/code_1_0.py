import math

# Part 1: Solve for the total number of pearls

# The number of pearls remaining on the string is "a seven shy of eleven times eleven"
pearls_on_string = (11 * 11) - 7

# The equation for the total number of pearls (x) is:
# x = x/6 + x/5 + x/3 + x/10 + pearls_on_string
# To solve for x, we rearrange the equation:
# x - (x/6 + x/5 + x/3 + x/10) = pearls_on_string
# Find a common denominator for the fractions (which is 30):
# x - (5x/30 + 6x/30 + 10x/30 + 3x/30) = pearls_on_string
# x - 24x/30 = pearls_on_string
# x - 4x/5 = pearls_on_string
# x/5 = pearls_on_string
# So, x = pearls_on_string * 5

total_pearls = pearls_on_string * 5

print("--- Part 1: Finding the total number of pearls ---")
print(f"First, let's verify the equation with the calculated total of {total_pearls} pearls:")

floor = total_pearls / 6
bed = total_pearls / 5
woman = total_pearls / 3
lover = total_pearls / 10

# Print each number in the final equation
print(f"The solved equation is: {math.ceil(floor)} + {math.ceil(bed)} + {math.ceil(woman)} + {math.ceil(lover)} + {pearls_on_string} = {total_pearls}")
print(f"So, the total number of pearls on the necklace was: {total_pearls}\n")


# Part 2: Solve for how many more pearls they need

print("--- Part 2: Finding how many more pearls are needed ---")

# The "fallen ones" are all pearls not on the string
fallen_pearls = total_pearls - pearls_on_string
print(f"Total pearls fallen (not on the string): {fallen_pearls}")

# They find 1/3rd of the fallen pearls
found_pearls = fallen_pearls / 3
print(f"Pearls they find back (1/3 of fallen): {math.ceil(found_pearls)}")

# The number they still need is the number of fallen pearls minus the ones they found
needed_pearls = fallen_pearls - found_pearls
print(f"Therefore, the number of pearls they still need is: {math.ceil(needed_pearls)}")

<<<304>>>