# Step 1: Define the topological properties of the base space N.
# N is a punctured torus, which is a torus with one disk removed.
# The genus of a punctured torus (g_N) is 1.
g_N = 1
# The number of boundary components of a punctured torus (k) is 1.
k = 1

# Step 2: Use the formula for the genus of a 2-sheeted branched cover M.
# The formula is g(M) = 2 * g(N) + k - 1.
g_M = 2 * g_N + k - 1

# Step 3: Print the calculation and the result.
print("The genus of the configuration space, g(M), is calculated using the formula:")
print(f"g(M) = 2 * g(N) + k - 1")
print(f"Substituting the values for a punctured torus (g(N)={g_N}, k={k}):")
print(f"g(M) = 2 * {g_N} + {k} - 1 = {g_M}")

# The final answer is the value of g_M
print("\nThe genus of the configuration space is:")
print(g_M)