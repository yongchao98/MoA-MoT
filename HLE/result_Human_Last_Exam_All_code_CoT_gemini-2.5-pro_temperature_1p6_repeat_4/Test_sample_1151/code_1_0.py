# Part 1: Calculate the total number of pearls on the necklace.

# The number of pearls remaining on the string is "a seven shy of eleven times eleven".
pearls_on_string = (11 * 11) - 7

# The total number of pearls, 'x', is the sum of its parts.
# The equation is: x = (1/6)x + (1/5)x + (1/3)x + (1/10)x + pearls_on_string
# To solve for x, we can rearrange the equation:
# x - (1/6 + 1/5 + 1/3 + 1/10)x = pearls_on_string
# The sum of the fractions of fallen pearls is:
fraction_fallen = (1/6) + (1/5) + (1/3) + (1/10)
# So, x * (1 - fraction_fallen) = pearls_on_string
# And, x = pearls_on_string / (1 - fraction_fallen)

total_pearls = pearls_on_string / (1 - fraction_fallen)
# The result should be a whole number, so we convert it to an integer.
total_pearls = int(total_pearls)

print("--- Part 1: How many pearls were there altogether? ---")
print("Let 'x' be the total number of pearls. The problem is described by the equation:")
print(f"x = (1/6)x + (1/5)x + (1/3)x + (1/10)x + {pearls_on_string}")
print(f"\nBy solving the equation, we find the total number of pearls on the necklace was: {total_pearls}\n")


# Part 2: Calculate how many more pearls they need.

# First, find the total number of pearls that fell off the string.
fallen_pearls = total_pearls - pearls_on_string

# Next, calculate the number of pearls they find back (1/3 of the fallen ones).
found_pearls = fallen_pearls / 3
found_pearls = int(found_pearls)

# The number of pearls they still need is the difference between what fell and what was found.
needed_pearls = fallen_pearls - found_pearls

print("--- Part 2: How many more pearls are needed? ---")
print(f"The total number of fallen pearls is: {total_pearls} - {pearls_on_string} = {fallen_pearls}")
print(f"They manage to find back 1/3 of the fallen pearls: {fallen_pearls} / 3 = {found_pearls}")
print(f"The number of pearls they still need to complete the necklace is: {fallen_pearls} - {found_pearls} = {needed_pearls}")