import math

# The final expression simplifies to (3/2) * 10^(10/3) + 37/4.
# Let's calculate the value of this expression.

coeff1 = 3/2
base = 10
exponent = 10/3
term2 = 37/4

# Calculate the first term
term1_val = coeff1 * (base ** exponent)

# Calculate the final result
result = term1_val + term2

# Print the final equation with the calculated numbers
# We show the components of the simplified expression
print(f"The calculation is based on the simplified expression: A * B^C + D")
print(f"Where A = {coeff1}, B = {base}, C = {exponent:.4f}, D = {term2}")
print(f"The equation with evaluated numbers is:")
print(f"{coeff1} * {base**exponent} + {term2} = {result}")
