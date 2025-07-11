import math

# Based on the analytical derivation, the integral simplifies to the expression:
# I = (8/15)*pi^8 + (1/3)*pi^2 - (1/2)*pi + 1

# Define pi
pi = math.pi

# The final equation is composed of four numerical terms.
# Let's define the numbers (coefficients and powers) in the final equation.
c1, p1 = 8/15, 8
c2, p2 = 1/3, 2
c3, p3 = -1/2, 1
c4 = 1

# Calculate the value of each term in the equation
term1 = c1 * (pi**p1)
term2 = c2 * (pi**p2)
term3 = c3 * (pi**p3)
term4 = c4

# The final result is the sum of these terms
result = term1 + term2 + term3 + term4

# As requested, here are the numbers in the final equation and the result.
print("The final simplified equation is:")
print(f"({c1:.4f}) * pi**{p1} + ({c2:.4f}) * pi**{p2} + ({c3:.2f}) * pi**{p3} + {c4}")
print("\nBreaking it down term by term:")
print(f"Term 1: (8/15) * pi^8 = {term1}")
print(f"Term 2: (1/3) * pi^2 = {term2}")
print(f"Term 3: (-1/2) * pi = {term3}")
print(f"Term 4: 1 = {term4}")
print("-------------------------------------------")
print(f"The final value of the integral is: {result}")