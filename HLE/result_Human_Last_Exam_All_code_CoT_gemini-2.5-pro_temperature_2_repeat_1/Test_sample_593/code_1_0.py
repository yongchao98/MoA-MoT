# Define the treewidths of graphs H and G, and the number of shared vertices.
# These are example values.
t_H = 5
t_G = 7
k = 4

# The tight upper bound for the treewidth of the combined graph F is max(t_H, t_G, k-1).
k_minus_1 = k - 1
treewidth_F_bound = max(t_H, t_G, k_minus_1)

# Print the final equation with the numbers plugged in.
# The 'tw' notation stands for treewidth.
print(f"Given t_H = {t_H}, t_G = {t_G}, k = {k}:")
print("The tight upper bound on the treewidth of F is given by the formula: tw(F) <= max(t_H, t_G, k-1)")
print(f"tw(F) <= max({t_H}, {t_G}, {k}-1)")
print(f"tw(F) <= max({t_H}, {t_G}, {k_minus_1})")
print(f"tw(F) <= {treewidth_F_bound}")