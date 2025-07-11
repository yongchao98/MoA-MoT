# Define the parameters based on the problem statement.

# N is the number of free groups in the free product, G = F_1 * ... * F_N.
N = 19
# M is the exponent of each commutator, c = product(c_i^M).
M = 30

# The stable commutator length (scl) of a single commutator c_i = [a_i, b_i]
# in the free group F_i is a known value of 1/2.
scl_ci = 0.5

# The total stable commutator length is calculated based on two key properties:
# 1. Additivity over free products: scl_G(product) = sum(scl_F_i).
# 2. Scaling with exponent: scl(g^M) = M * scl(g).
# Combining these, the total scl is: N * M * scl_ci.
total_scl = N * M * scl_ci

# Format the result as an integer if it's a whole number.
if total_scl.is_integer():
    total_scl = int(total_scl)

# Print the final equation, showing each number used in the calculation.
print(f"The final calculation for the stable commutator length is:")
print(f"{N} * {M} * {scl_ci} = {total_scl}")