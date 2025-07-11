# Parameters for the graphs H and G
# t_H: treewidth of H
# t_G: treewidth of G
# k: number of vertices in the intersection of V(H) and V(G)

# We use the values from the example demonstrating the tightness of the bound.
# H is a 5-cycle, G is a 5-cycle, their union is a K_5.
t_H = 2
t_G = 2
k = 5

# The tight upper bound on the treewidth of the union graph F is max(t_H, t_G, k - 1).
# Here, we calculate this bound.
result = max(t_H, t_G, k - 1)

# We print the final equation showing the calculation.
print("A tight upper bound for the treewidth of F = H U G is given by max(t_H, t_G, k - 1).")
print("Using the example values where H and G are 5-cycles whose union is K_5:")
print(f"t_H = {t_H}")
print(f"t_G = {t_G}")
print(f"k = {k}")
print(f"The calculation is: max({t_H}, {t_G}, {k} - 1) = {result}")
