import numpy as np

# The group G is the direct product Z_4 x Z_2 x Z_2.
# Elements are represented as numpy arrays [a, b, c] where a is in {0,1,2,3} and b,c are in {0,1}.
group_orders = np.array([4, 2, 2])

def add(g1, g2):
    """Adds two group elements."""
    return (np.array(g1) + np.array(g2)) % group_orders

def double(g):
    """Doubles a group element."""
    return (2 * np.array(g)) % group_orders

# A maximal sum-free set S in G = Z_4 x Z_2 x Z_2
S = [[2, 0, 0], [0, 1, 0], [2, 1, 1]]

# Calculate the size of S
size_S = len(S)

# Calculate the set k(S)
k_S = []
# Iterate through all 16 elements of the group G
for a in range(4):
    for b in range(2):
        for c in range(2):
            g = [a, b, c]
            doubled_g = double(g)
            
            # Check if 2g is in S
            for s_elem in S:
                if np.array_equal(doubled_g, s_elem):
                    k_S.append(g)
                    break # Avoid adding g multiple times if it maps to different s

# Calculate the size of k(S)
size_k_S = len(k_S)

# The problem is to find the smallest |G| with a maximal sum-free S where |k(S)| > 2|S|.
# Our example uses a group of size 16.
group_size = np.prod(group_orders)

print(f"The group G = Z_4 x Z_2 x Z_2 has size |G| = {group_size}.")
print(f"A maximal sum-free set is S = {S}.")
print(f"The size of S is |S| = {size_S}.")
print(f"The set k(S) = {{g in G | 2g in S}} is {k_S}.")
print(f"The size of k(S) is |k(S)| = {size_k_S}.")
print("\nChecking the condition |k(S)| > 2|S|:")
print(f"{size_k_S} > 2 * {size_S}")
print(f"{size_k_S} > {2 * size_S}")
print(f"The condition is {'met' if size_k_S > 2 * size_S else 'not met'}.")

print(f"\nThe smallest size of such a group is {group_size}.")

print("\nTo show the final numbers in the equation: ")
print(f"{size_k_S} > 2 * {size_S}")

<<<16>>>