import math

# Step 1: Define the optimal values for the parameters derived from the analysis.
# We found that the maximum of |b| + |c| is achieved when |b| and |c| are
# related in a specific way. Let's use the positive values for b and c that
# lead to the maximum.
# k0 corresponds to |b|
k0 = 2 / math.sqrt(5)
# u corresponds to |c|
u = (math.sqrt(5) + 1) / (2 * math.sqrt(5))

# These values correspond to coefficients of a valid polynomial f(x) = ax^2+bx+c.
# For example, one such polynomial has:
# b_opt = k0
# c_opt = u
# a_opt = -u
# For this polynomial, the sum |b| + |c| is k0 + u.
b_val = k0
c_val = u

# Step 2: Calculate the maximum value of |b| + |c|
max_value = abs(b_val) + abs(c_val)

# Step 3: Output the numbers in the final equation as requested.
# The equation is |b| + |c| = max_value
print("An optimal choice for |b| and |c| gives the equation:")
print(f"{abs(b_val)} + {abs(c_val)} = {max_value}")
print("\nThe maximum value is the golden ratio, phi.")
print(f"Maximum value: {max_value}")

# For verification, the value is (sqrt(5)+1)/2
phi = (math.sqrt(5) + 1) / 2
print(f"Value of phi for comparison: {phi}")