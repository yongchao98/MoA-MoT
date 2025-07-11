# Step 1: Define the dimension of the space.
# The curve is in R^3, so the dimension d is 3.
d = 3

# Step 2: Define the components of the formula.
# According to the Bourgain-Demeter-Guth theorem, the sharp l^2 decoupling
# exponent for a non-degenerate curve in R^d is p = d * (d + 1).
term1 = d
term2 = d + 1

# Step 3: Calculate the exponent.
exponent = term1 * term2

# Step 4: Print the final equation with all its components and the result.
print(f"The sharp l^2 decoupling exponent p is calculated using the formula: p = d * (d + 1)")
print(f"For the given curve, the dimension d is {d}.")
print(f"The calculation is: {term1} * ({term1} + 1) = {exponent}")