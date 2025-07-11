# Let's define example values for the treewidth of graph H (t_H),
# the treewidth of graph G (t_G), and the number of shared vertices (k).

# The treewidth of graph H.
t_H = 8

# The treewidth of graph G.
t_G = 6

# The number of vertices in the intersection of H and G.
k = 10

# The tight upper bound for the treewidth of the combined graph F is
# given by the formula: max(t_H, t_G, k - 1).

# First, calculate the third term in the max function.
k_minus_1 = k - 1

# Now, calculate the final upper bound.
upper_bound = max(t_H, t_G, k_minus_1)

# Finally, print the result in a clear format showing the full equation.
print("--- Calculating the Treewidth Upper Bound for a Graph Union ---")
print(f"Given the treewidth of H, t_H = {t_H}")
print(f"Given the treewidth of G, t_G = {t_G}")
print(f"Given the number of shared vertices, k = {k}")
print("\nThe tight upper bound is calculated as: max(t_H, t_G, k - 1)")
print("\nSubstituting the given values into the equation:")
print(f"max({t_H}, {t_G}, {k} - 1) = max({t_H}, {t_G}, {k_minus_1}) = {upper_bound}")
print("\nTherefore, the tight upper bound on the treewidth of F is", upper_bound)
