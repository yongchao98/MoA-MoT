# Set the dimensions of the two spheres from the problem.
n1 = 4
n2 = 6

# The connectivity of an n-sphere S^n is (n-1).
# This is because the first non-trivial homotopy group of S^n is π_n(S^n).
# So π_i(S^n) is zero for i < n.
conn_s_n1 = n1 - 1
conn_s_n2 = n2 - 1

# According to a theorem in homotopy theory, the connectivity of the map
# Σ(ΩX ∧ ΩY) → Ω(X ∧ Y) for suspension spaces X and Y is given by
# conn(X) + conn(Y) + 1.
# Since S^4 and S^6 are suspensions (S^4 = ΣS^3, S^6 = ΣS^5), we can apply this formula.
map_connectivity = conn_s_n1 + conn_s_n2 + 1

# Print the final equation with all its components.
print(f"{conn_s_n1} + {conn_s_n2} + 1 = {map_connectivity}")
