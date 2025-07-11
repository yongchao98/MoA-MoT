import math

def n_cycles_in_K(n):
    """Calculates the total number of simple cycles in a complete graph K_n."""
    if n < 3:
        return 0
    total_cycles = 0
    for k in range(3, n + 1):
        # Number of ways to choose k vertices: C(n, k)
        # Number of ways to form a cycle with k vertices: (k-1)!/2
        cycles_k = math.comb(n, k) * math.factorial(k - 1) // 2
        total_cycles += cycles_k
    return total_cycles

def n_cycles_with_edge_in_K(n):
    """Calculates the number of cycles containing a specific edge in K_n."""
    if n < 3:
        return 0
    total_cycles_with_edge = 0
    for k in range(3, n + 1):
        # To form a k-cycle with a given edge (u,v), we need to choose k-2 other vertices
        # from the remaining n-2 vertices and form a path of length k-1 between u and v.
        # The number of such paths is P(n-2, k-2).
        if n - 2 >= k - 2:
            cycles_k_with_edge = math.perm(n - 2, k - 2)
            total_cycles_with_edge += cycles_k_with_edge
    return total_cycles_with_edge
    
def n_paths_between_endpoints_in_K_minus_e(n):
    """
    Calculates the number of simple paths between two vertices u, v in K_n,
    assuming the direct edge (u,v) has been removed.
    """
    if n < 2:
        return 0
    total_paths = 0
    # Paths can have length l from 2 to n-1.
    for l in range(2, n):
        # A path of length l has l-1 intermediate vertices.
        # We choose l-1 vertices from the n-2 available vertices and order them.
        if n-2 >= l-1:
            paths_l = math.perm(n - 2, l - 1)
            total_paths += paths_l
    return total_paths

# --- Step 1: Calculate cycles in the base cliques K_5 and K_4 ---
cycles_k5 = n_cycles_in_K(5)
cycles_k4 = n_cycles_in_K(4)
print(f"A complete graph K_5 has {cycles_k5} cycles.")
print(f"A complete graph K_4 has {cycles_k4} cycles.")
print("-" * 20)

# --- Step 2: Calculate cycles lost by removing one edge from each clique ---
cycles_lost_k5 = n_cycles_with_edge_in_K(5)
cycles_lost_k4 = n_cycles_with_edge_in_K(4)
cycles_in_k5_minus_e = cycles_k5 - cycles_lost_k5
cycles_in_k4_minus_e = cycles_k4 - cycles_lost_k4

print(f"Removing an edge from K_5 destroys {cycles_lost_k5} cycles.")
print(f"The remaining K_5-minus-edge component has {cycles_in_k5_minus_e} cycles.")
print(f"Removing an edge from K_4 destroys {cycles_lost_k4} cycles.")
print(f"The remaining K_4-minus-edge component has {cycles_in_k4_minus_e} cycles.")
print("-" * 20)

# --- Step 3: Calculate the number of new cycles formed ---
# New cycles are formed by: Path in (K_5-e) -> new_edge_1 -> Path in (K_4-e) -> new_edge_2 -> back
# The number of new cycles is the product of the number of paths between the endpoints of the new edges.

# Paths in the K_5 component between the endpoints of the removed edge
paths_in_k5_comp = n_paths_between_endpoints_in_K_minus_e(5)

# Paths in the K_4 component. This is a special case calculation.
# The graph is K_4 minus an edge, say (v3,v4). We need the number of paths
# between two other vertices, say v1 and v2.
# Path v1-v2 (direct edge): 1
# Path v1-v3-v2 (length 2): 1
# Path v1-v4-v2 (length 2): 1
# So, total paths = 1 + 1 + 1 = 3.
# Let's assume the connecting edges are (u1,v1) and (u2,v2) where u1,u2 are from K5-e and v1,v2 from K4-e
# and e1=(u1,u2) and e2=(v1,v2). We need paths from u1 to u2 and v1 to v2.
paths_in_k4_comp = 3 
# Note: For the specific construction that maximizes cycles, the number of paths in the K4-e component is 3.

new_cycles = paths_in_k5_comp * paths_in_k4_comp

print(f"The number of paths between the connection points in the K_5 component is {paths_in_k5_comp}.")
print(f"The number of paths between the connection points in the K_4 component is {paths_in_k4_comp}.")
print(f"The number of new cycles formed by connecting the components is {paths_in_k5_comp} * {paths_in_k4_comp} = {new_cycles}.")
print("-" * 20)

# --- Step 4: Calculate the total maximum number of cycles ---
total_max_cycles = cycles_in_k5_minus_e + cycles_in_k4_minus_e + new_cycles

print("The final equation for the total number of Mexican standoffs is:")
print(f"Total = (Cycles in K5 - Cycles lost) + (Cycles in K4 - Cycles lost) + (New Cycles)")
print(f"Total = ({cycles_k5} - {cycles_lost_k5}) + ({cycles_k4} - {cycles_lost_k4}) + ({paths_in_k5_comp} * {paths_in_k4_comp})")
print(f"Total = {cycles_in_k5_minus_e} + {cycles_in_k4_minus_e} + {new_cycles}")
print(f"The maximum number of standoffs is: {total_max_cycles}")