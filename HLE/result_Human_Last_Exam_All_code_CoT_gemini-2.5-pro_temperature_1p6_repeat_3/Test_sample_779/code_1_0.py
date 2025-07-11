import collections

def check_g_admissible_coarsening(graph_n, edges, part_rho, part_sigma):
    """
    Checks if partition sigma is a G-admissible coarsening of partition rho.
    
    A G-admissible coarsening is defined as being obtained by merging exactly two
    distinct blocks B_i and B_j of rho, where there exist vertices v_i in B_i
    and v_j in B_j with {v_i, v_j} being an edge in the graph G.
    """
    
    # Sigma must have exactly one less block than rho.
    if len(part_sigma) != len(part_rho) - 1:
        return False
        
    # Convert partitions to a canonical form (frozenset of frozensets of ints) for easy comparison.
    rho_fset = frozenset(frozenset(b) for b in part_rho)
    sigma_fset = frozenset(frozenset(b) for b in part_sigma)

    # Find the single block in sigma that is not in rho. This must be the merged block.
    merged_block = None
    sigma_only_blocks = sigma_fset.difference(rho_fset)
    if len(sigma_only_blocks) != 1:
        return False
    merged_block = set(list(sigma_only_blocks)[0])

    # Find the two blocks in rho that are not in sigma. These must be the original blocks.
    rho_only_blocks = rho_fset.difference(sigma_fset)
    if len(rho_only_blocks) != 2:
        return False
    
    original_blocks = [set(b) for b in rho_only_blocks]
    b_i, b_j = original_blocks[0], original_blocks[1]
    
    # Verify that merging b_i and b_j creates the new block in sigma.
    if b_i.union(b_j) != merged_block:
        return False

    # Check for G-admissibility: there must be an edge in G connecting b_i and b_j.
    for v_i in b_i:
        for v_j in b_j:
            # Check for edge in both directions since the edge set might not be canonical.
            if (v_i, v_j) in edges or (v_j, v_i) in edges:
                # Found a connecting edge, the merge is G-admissible.
                return True
    
    # If no connecting edge is found after checking all pairs, it's not G-admissible.
    return False

# ---- Example Usage ----
# Let G be the cycle graph C4 on 4 vertices. This graph is not chordal.
n = 4
# Edges are represented as a set of tuples for efficient lookup.
E = {(1, 2), (2, 3), (3, 4), (4, 1)}

# Let rho be the all-singleton partition, bot_4
rho = [{1}, {2}, {3}, {4}]
# Let sigma be a partition where {1} and {2} are merged.
sigma = [{1, 2}, {3}, {4}]

# This merge is admissible because the edge {1, 2} exists in C4.
is_admissible = check_g_admissible_coarsening(n, E, rho, sigma)

print(f"Let G be the cycle graph C4 with V(G) = {{1, 2, 3, 4}}.")
print(f"Let rho = {rho}")
print(f"Let sigma = {sigma}")
print(f"Is sigma a G-admissible coarsening of rho? {is_admissible}")

# Another example: try to merge {1} and {3}.
sigma_2 = [{1, 3}, {2}, {4}]
# This merge is not admissible because there is no edge {1, 3} in C4.
is_admissible_2 = check_g_admissible_coarsening(n, E, rho, sigma_2)
print(f"\nLet sigma_2 = {sigma_2}")
print(f"Is sigma_2 a G-admissible coarsening of rho? {is_admissible_2}")
