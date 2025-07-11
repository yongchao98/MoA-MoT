import math

# Define the value of n
n = 1000

# Constants in the final equation
const_4 = 4
const_2 = 2

# Calculate the denominator in the argument of cosine
denominator = n + const_2

# Calculate the final value based on the derived formula 1/p_n = 4 * cos^2(pi / (n + 2))
result = const_4 * (math.cos(math.pi / denominator) ** 2)

# Print the components of the final equation and the result
print("The final equation is 1/p_1000 = 4 * cos^2(pi / (1000 + 2))")
print("The number", const_4, "multiplied by the square of cosine of (pi divided by", denominator, ") equals the result.")
print("Result:", result)