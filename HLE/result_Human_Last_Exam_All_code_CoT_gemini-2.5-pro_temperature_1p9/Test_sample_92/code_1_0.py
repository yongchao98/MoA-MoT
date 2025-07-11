# Define the parameters of the problem
escape_bin = 2025
melt_bin = 2024
start_bin = 0

# For a symmetric 1D random walk, the probability of absorption at one boundary
# is a linear function of the starting position.
# Let p_n be the probability of escaping starting from bin n.
# p_n = A*n + B
# Using the boundary conditions:
# p_escape = A * escape_bin + B = 1
# p_melt = A * melt_bin + B = 0

# Solving for A and B:
# A = (1 - 0) / (escape_bin - melt_bin)
# B = -A * melt_bin

A = 1 / (escape_bin - melt_bin)
B = -A * melt_bin

# The probability of escaping from the start_bin is p_start:
p_start = A * start_bin + B

# We can also write this as:
# p_start = (start_bin - melt_bin) / (escape_bin - melt_bin)
numerator = start_bin - melt_bin
denominator = escape_bin - melt_bin
result = numerator / denominator

print(f"The starting bin is {start_bin}.")
print(f"The escape bin is {escape_bin}.")
print(f"The melting bin is {melt_bin}.")
print("\nThe probability is calculated as (start_bin - melt_bin) / (escape_bin - melt_bin).")
print(f"So the equation is ({start_bin} - {melt_bin}) / ({escape_bin} - {melt_bin}).")
print(f"Result: {result}")