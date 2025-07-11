import networkx as nx
import math

def solve_subgraph_counting():
    """
    This function calculates the number of subgraphs of the Gosset graph
    that are isomorphic to the HoG graph ID 50698 (Kneser graph K(8,2)).
    """
    # 1. Construct the two graphs.
    # G1 is the larger graph: the Gosset graph.
    G1 = nx.gosset_graph()

    # G2 is the pattern graph: HoG graph ID 50698.
    # The problem description identifies it as the Kneser graph K(n,k)
    # where n=8 (from K8) and k=2 (edges of K8).
    n_kneser = 8
    k_kneser = 2
    G2 = nx.kneser_graph(n_kneser, k_kneser)

    print(f"Searching for subgraphs of the Gosset graph ({G1.number_of_nodes()} vertices, {G1.number_of_edges()} edges)")
    print(f"isomorphic to the Kneser graph K({n_kneser},{k_kneser}) ({G2.number_of_nodes()} vertices, {G2.number_of_edges()} edges).\n")
    
    # 2. Count the total number of isomorphic mappings from G2 to subgraphs of G1.
    # This process can be slow, so we'll print a status message.
    print("Counting all possible isomorphic mappings... (this may take a moment)")
    graph_matcher = nx.isomorphism.GraphMatcher(G1, G2)
    num_isomorphisms = sum(1 for _ in graph_matcher.subgraph_isomorphisms_iter())
    
    # 3. Calculate the size of the automorphism group of the pattern graph G2.
    # For a Kneser graph K(n, k) with n > 2k, the size of its automorphism
    # group is n!. Here, 8 > 2*2=4, so the formula applies.
    automorphism_group_size = math.factorial(n_kneser)

    # 4. Calculate the number of unique subgraphs.
    # This is the total number of mappings divided by the number of symmetries
    # (automorphisms) of the pattern graph. The result must be an integer.
    if num_isomorphisms == 0:
        num_subgraphs = 0
    else:
        if num_isomorphisms % automorphism_group_size != 0:
            print("Error: The number of isomorphisms is not divisible by the automorphism group size.")
            # Fallback to float division if something is unexpected
            num_subgraphs = num_isomorphisms / automorphism_group_size
        else:
            # Use integer division for a clean result
            num_subgraphs = num_isomorphisms // automorphism_group_size

    # 5. Print the final calculation as requested.
    print("\nFinal Calculation:")
    print("------------------")
    print(f"Total Isomorphic Mappings Found: {num_isomorphisms}")
    print(f"Symmetries of Pattern Graph (size of Aut(K(8,2)) = 8!): {automorphism_group_size}")
    print("\nNumber of unique subgraphs = (Total Isomorphic Mappings) / (Symmetries of Pattern Graph)")
    print(f"Result = {num_isomorphisms} / {automorphism_group_size} = {num_subgraphs}")

if __name__ == '__main__':
    solve_subgraph_counting()
