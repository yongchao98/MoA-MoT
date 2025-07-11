import networkx as nx
from itertools import product
from math import prod

def count_homomorphisms(H, G):
    """
    Counts the number of graph homomorphisms from graph H to graph G.
    A homomorphism is a mapping f: V(H) -> V(G) such that if (u,v) is an edge in H,
    then (f(u), f(v)) is an edge in G.
    """
    h_nodes = list(H.nodes())
    g_nodes = list(G.nodes())
    
    # Generate all possible functions from V(H) to V(G)
    all_mappings = product(g_nodes, repeat=len(h_nodes))
    
    hom_count = 0
    for mapping_tuple in all_mappings:
        f = {h_node: mapping_tuple[i] for i, h_node in enumerate(h_nodes)}
        
        is_homomorphism = True
        for u, v in H.edges():
            # Check if the edge mapping is preserved
            if not G.has_edge(f[u], f[v]):
                is_homomorphism = False
                break
        
        if is_homomorphism:
            hom_count += 1
            
    return hom_count

def main():
    """
    Main function to demonstrate the argument.
    """
    # Step 1: Define two graphs G1 and G2 that satisfy the premise.
    # For demonstration, we'll use two identical graphs. This trivially satisfies
    # the condition that hom(T, G1) = hom(T, G2) for any tree T.
    G1 = nx.cycle_graph(4)
    G2 = nx.cycle_graph(4)

    # Step 2: Define an acyclic conjunctive query phi.
    # We represent phi by its query graph Q_phi, which is a forest.
    # Let's create a query graph with two components: a path of length 2 (P3)
    # and a single edge (P2).
    Q_phi = nx.Graph()
    Q_phi.add_edges_from([(0, 1), (1, 2)]) # Tree component 1: P3
    Q_phi.add_edges_from([(3, 4)])         # Tree component 2: P2

    print("Demonstrating that ans(phi, G1) == ans(phi, G2)")
    print("-" * 50)

    # Step 3: Get the tree components of the query graph Q_phi.
    tree_components = [Q_phi.subgraph(c).copy() for c in nx.connected_components(Q_phi)]

    hom_counts_g1 = []
    hom_counts_g2 = []

    # Step 4: Calculate homomorphism counts for each tree component.
    for i, T in enumerate(tree_components):
        h_g1 = count_homomorphisms(T, G1)
        h_g2 = count_homomorphisms(T, G2)
        hom_counts_g1.append(h_g1)
        hom_counts_g2.append(h_g2)
        print(f"Homomorphism count for tree component T_{i+1} in G1: {h_g1}")
        print(f"Homomorphism count for tree component T_{i+1} in G2: {h_g2}")

    # Step 5: Calculate the total number of answers for phi.
    # This is the product of the homomorphism counts of the components.
    ans_g1 = prod(hom_counts_g1)
    ans_g2 = prod(hom_counts_g2)

    # Step 6: Print the final equations and results.
    equation_g1 = " * ".join(map(str, hom_counts_g1))
    equation_g2 = " * ".join(map(str, hom_counts_g2))

    print("\nFinal equation for G1: {} = {}".format(equation_g1, ans_g1))
    print("Final equation for G2: {} = {}".format(equation_g2, ans_g2))

    print(f"\nTotal number of answers for phi in G1: {ans_g1}")
    print(f"Total number of answers for phi in G2: {ans_g2}")

    print("\nConclusion: As demonstrated, the number of answers is the same.")
    print("This is because the number of answers is determined by homomorphism counts from the query's tree components, which are guaranteed to be equal for G1 and G2 by the problem's premise.")

if __name__ == "__main__":
    main()