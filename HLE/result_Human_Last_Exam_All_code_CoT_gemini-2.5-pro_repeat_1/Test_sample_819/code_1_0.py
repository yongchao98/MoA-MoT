import networkx as nx

def get_graph_properties(G):
    """Helper function to get degrees and number of edges."""
    degrees = [val for (node, val) in G.degree()]
    num_edges = G.number_of_edges()
    return degrees, num_edges

def hom_P2(G):
    """Computes the number of homomorphisms from P_2 (a single edge) to G."""
    # A homomorphism from an edge (u,v) to G corresponds to picking an edge
    # in G and mapping (u,v) to it in one of two ways.
    # So, hom(P_2, G) = 2 * |E(G)|.
    return 2 * G.number_of_edges()

def hom_P3(G):
    """Computes the number of homomorphisms from P_3 (path of length 2) to G."""
    # A homomorphism from P_3 (u-v-w) to G is a walk of length 2.
    # The number of such walks is sum(deg(c)^2 for c in V(G)).
    return sum(d**2 for d in dict(G.degree()).values())

def hom_P4(G):
    """Computes the number of homomorphisms from P_4 (path of length 3) to G."""
    # A homomorphism from P_4 is a walk of length 3.
    # The number is Sum_{u,v} (A^3)_{uv} = 1^T A^3 1.
    # This can be calculated by iterating through all walks.
    count = 0
    for u in G.nodes():
        for v in G.neighbors(u):
            for w in G.neighbors(v):
                for x in G.neighbors(w):
                    count += 1
    return count

def main():
    """
    Main function to demonstrate the conclusion.
    """
    # 1. Define two simple graphs, G1 and G2.
    # For the premise hom(T,G1) = hom(T,G2) to hold for all trees T,
    # G1 and G2 must be isomorphic.
    # Let's use two isomorphic 5-cycles (C5).
    G1 = nx.cycle_graph(5)  # Nodes 0, 1, 2, 3, 4
    
    # G2 is isomorphic to G1, just with different node labels.
    G2 = nx.Graph()
    G2.add_edges_from([('a','b'), ('b','c'), ('c','d'), ('d','e'), ('e','a')])

    # The premise holds because G1 and G2 are isomorphic.
    # We can verify this for a few trees.
    hom_p2_g1 = hom_P2(G1)
    hom_p2_g2 = hom_P2(G2)
    hom_p3_g1 = hom_P3(G1)
    hom_p3_g2 = hom_P3(G2)

    print(f"Homomorphisms from P_2 to G1: {hom_p2_g1}")
    print(f"Homomorphisms from P_2 to G2: {hom_p2_g2}")
    print(f"Homomorphisms from P_3 to G1: {hom_p3_g1}")
    print(f"Homomorphisms from P_3 to G2: {hom_p3_g2}")
    print("-" * 20)

    # 2. Define an acyclic conjunctive query phi.
    # Let phi(x,y,z,w) be E(x,y) AND E(y,z) AND E(z,w).
    # This is a "full" query (no existential quantifiers).
    # Its query graph Q is the path graph P_4, which is a tree.
    # The number of answers of phi in a graph G is hom(P_4, G).
    
    # 3. Calculate the number of answers for G1 and G2.
    num_answers_g1 = hom_P4(G1)
    num_answers_g2 = hom_P4(G2)

    print("Query phi(x,y,z,w) = E(x,y) AND E(y,z) AND E(z,w)")
    print("The query graph Q is a path P_4.")
    print(f"Number of answers of phi in G1 = hom(P_4, G1) = {num_answers_g1}")
    print(f"Number of answers of phi in G2 = hom(P_4, G2) = {num_answers_g2}")
    print("-" * 20)

    # 4. State the conclusion.
    print("Conclusion:")
    print("The number of answers is the same for both graphs.")
    print("This is because the premise hom(T, G1) = hom(T, G2) for all trees T implies that G1 and G2 must be isomorphic.")
    print("Isomorphic graphs have the same number of answers for any query.")
    print("Therefore, it is not possible for G1 and G2 to have different numbers of answers.")


if __name__ == "__main__":
    main()
