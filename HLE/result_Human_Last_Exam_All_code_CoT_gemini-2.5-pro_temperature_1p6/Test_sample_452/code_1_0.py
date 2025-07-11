import math

# The problem is to find the constant 'b' in the asymptotic formula
# C(n) ~ b * n^(3/2) for the expected cover-and-return time on a random tree.

# Step 1: The expected cover-and-return time, C(T), averaged over starting vertices,
# is related to the tree's Wiener index W(T).
# For large n, C(T) ~ (4/n) * W(T).
# We need the expectation over random trees, so C(n) ~ (4/n) * E[W(T_n)].

# Step 2: A known result for the expected Wiener index of a random n-vertex tree is:
# E[W(T_n)] ~ (sqrt(2*pi)/8) * n^(5/2)
# So we define the coefficient for the Wiener index formula.
E_W_coeff = math.sqrt(2 * math.pi) / 8

# Step 3: We substitute this into the expression for C(n).
# C(n) ~ (4/n) * (sqrt(2*pi)/8) * n^(5/2)
# C(n) ~ 4 * (sqrt(2*pi)/8) * n^(3/2)
# C(n) ~ (sqrt(2*pi)/2) * n^(3/2)

# Step 4: The constant b is therefore sqrt(2*pi)/2, which simplifies to sqrt(pi/2).
b = math.sqrt(math.pi / 2)

# The final result is the equation with the calculated value of b.
# Let's print out the key components.
print(f"The asymptotic formula for the expected Wiener Index is E[W(n)] ~ k * n^(5/2)")
print(f"where k = sqrt(2*pi)/8 = {E_W_coeff}")
print(f"The expected cover-and-return time C(n) is asymptotically related by C(n) ~ (4/n) * E[W(n)]")
print(f"Therefore, C(n) ~ (4/n) * ({E_W_coeff}) * n^(5/2)")
print(f"C(n) ~ {4*E_W_coeff} * n^(3/2)")
print(f"The constant b is sqrt(pi/2)")
print(f"The final equation is C(n) ~ {b} * n^(3/2)")

<<<1.2533141373155003>>>