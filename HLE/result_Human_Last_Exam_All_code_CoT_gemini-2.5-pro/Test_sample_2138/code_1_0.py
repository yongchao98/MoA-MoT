import math

# The integral J is found to be the product of two terms derived from the theory of special functions.
# The first term comes from the Gamma function part of the evaluation.
term1 = 2 * math.sqrt(math.pi)

# The second term is the value of the Hurwitz zeta function part, given by a known identity.
term2 = 2 * math.sqrt(2 * math.pi / 3)

# The final result is the product of these two terms.
result = term1 * term2

# We present the final equation with its components.
print(f"The analytical value is the result of the following product:")
# Using string formatting to show the equation with the numerical values.
# Note: the numbers are displayed with a certain precision for readability.
print(f"{term1:.6f} * {term2:.6f} = {result:.6f}")
print("The exact analytical value is 4 * pi * sqrt(6) / 3.")
# Final equation to be printed as requested by the user prompt
print(f"The final equation is: 2 * sqrt(pi) * 2 * sqrt(2*pi/3) = 4*pi*sqrt(6)/3")
print(f"Value = {result}")
