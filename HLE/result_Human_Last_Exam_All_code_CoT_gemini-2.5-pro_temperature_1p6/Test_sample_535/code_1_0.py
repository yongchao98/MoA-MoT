import math

# Define the constants from the characteristic equation r^2 - r - 1 = 0
sqrt5 = math.sqrt(5)
r1 = (1 + sqrt5) / 2
r2 = (1 - sqrt5) / 2

# The expression to be calculated is (1/sqrt(5)) * ( (10^10)^r1 + (10^10)^r2 )
# which simplifies to (1/sqrt(5)) * ( 10^(10*r1) + 10^(10*r2) )
exp1 = 10 * r1
exp2 = 10 * r2
term1 = 10**exp1
term2 = 10**exp2

# Calculate the final value
value = (1 / sqrt5) * (term1 + term2)

# Print the final equation with all its components
print(f"Result = (1/√5) * (10^(10 * ({r1:.5f})) + 10^(10 * ({r2:.5f})))")
print(f"Result = (1/{sqrt5:.5f}) * (10^{exp1:.5f} + 10^{exp2:.5f})")
print(f"Result = (1/{sqrt5:.5f}) * ({term1:.5e} + {term2:.5e})")
print(f"Result = {value}")
