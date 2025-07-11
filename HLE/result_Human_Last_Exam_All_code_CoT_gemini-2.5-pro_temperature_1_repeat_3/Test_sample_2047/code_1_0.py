import math

# Define parameters d and lambda.
# As per the problem description, d >= 4 and lambda >= 1.
# We will use d=4 and lambda=1.0 for this example.
d = 4
lmbda = 1.0

# The derived formula for l(d, lambda) is:
# l(d, lambda) = (1 / (2 * lambda)) * [arccos(sqrt(2/d))^2 - arccos(sqrt(3/d))^2]

# Calculate the terms inside the arccos functions
term_x1 = 3.0 / d
term_x2 = 2.0 / d

# Calculate the squared angles (geodesic distances squared)
w1_sq = math.acos(math.sqrt(term_x1))**2
w2_sq = math.acos(math.sqrt(term_x2))**2

# Calculate the final result
result = (1.0 / (2.0 * lmbda)) * (w2_sq - w1_sq)

# Print the final equation with all intermediate values
print("The final equation is:")
print(f"l({d}, {lmbda}) = (1 / (2 * {lmbda})) * [arccos(sqrt(2/{d}))^2 - arccos(sqrt(3/{d}))^2]")
print("Substituting values:")
print(f"l({d}, {lmbda}) = (1 / ({2 * lmbda})) * [arccos({math.sqrt(term_x2):.4f})^2 - arccos({math.sqrt(term_x1):.4f})^2]")
print(f"l({d}, {lmbda}) = (1 / ({2 * lmbda})) * [{math.acos(math.sqrt(term_x2)):.4f}^2 - {math.acos(math.sqrt(term_x1)):.4f}^2]")
print(f"l({d}, {lmbda}) = ({1.0/(2*lmbda):.4f}) * [{w2_sq:.4f} - {w1_sq:.4f}]")
print(f"l({d}, {lmbda}) = ({1.0/(2*lmbda):.4f}) * [{(w2_sq - w1_sq):.4f}]")
print(f"l({d}, {lmbda}) = {result:.4f}")

# The final answer in the required format
# <<<result>>>