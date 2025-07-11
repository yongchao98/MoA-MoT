import math

# The problem reduces to finding the integer part of 10^4 * L,
# where L is the limit L = lim_{N->inf} F(N)/ln(N).
# Our derivation shows that L = 2 / ln(alpha), where alpha is a root of r^2 - 13r + 1 = 0.

# The numbers in the final equation for the limit L = A / ln((B + sqrt(C))/D)
A = 2
B = 13
C = 165
D = 2

# We print the numbers as requested.
print(f"The limit L is determined by the equation: L = {A} / ln(({B} + sqrt({C}))/{D})")
print(f"The numbers in the equation are: A={A}, B={B}, C={C}, D={D}")

# Calculate alpha
alpha = (B + math.sqrt(C)) / D

# Calculate the limit L
limit_L = A / math.log(alpha)

# The value requested in the problem
result = 10**4 * limit_L

# Print the final integer part
print(f"\nThe calculated value of 10^4 * L is: {result}")
print(f"The integer part is: {math.floor(result)}")
