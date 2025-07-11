import math

# Define the value of n
n = 1000

# The formula for p_n is 1 / (4 * cos^2(pi / (n + 2)))
# We need to find the value of 1 / p_n for n = 1000.
# So we need to calculate 4 * cos^2(pi / (1000 + 2))

# The numbers in the final equation
four = 4
denominator = n + 2

# Calculate the result
result = four * (math.cos(math.pi / denominator))**2

# Print the final equation with the computed value
print(f"{four} * (cos(pi / {denominator}))^2 = {result}")