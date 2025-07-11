import networkx as nx
import itertools

def demonstrate_maximal_planar_properties():
    """
    This function uses a concrete example to illustrate the properties of adding
    an edge to a maximal planar graph.
    """
    
    # Let's construct G, a maximal planar graph on n=5 vertices.
    # A maximal planar graph on n=5 vertices has 3*n - 6 = 9 edges.
    # The complete graph K5 has 10 edges and is not planar.
    # A graph formed by removing one edge from K5 is maximal planar.
    # We will use vertices 0, 1, 2, 3, 4.
    G = nx.complete_graph(5)
    e_to_remove = (3, 4)
    G.remove_edge(*e_to_remove)
    
    # The edge 'e' that is not in G is the one we removed.
    e_to_add = e_to_remove
    
    # G' is G union {e}, which brings us back to K5.
    G_prime = G.copy()
    G_prime.add_edge(*e_to_add)

    print("Step 1: Setup")
    print(f"Let G be the graph K5 - e, where e = {e_to_remove}. This is a maximal planar graph.")
    print(f"Let G' = G + e, which is K5.")
    print("-" * 30)

    print("Step 2: Verify planarity properties")
    is_G_planar = nx.is_planar(G)
    print(f"Is G planar? {is_G_planar}")
    
    is_G_prime_planar = nx.is_planar(G_prime)
    print(f"Is G' = G + e planar? {is_G_prime_planar}")
    print("As expected, adding the edge 'e' to the maximal planar graph G makes it non-planar.")
    print("This implies any drawing of G' must have at least one crossing.")
    print("-" * 30)
    
    print("Step 3: Analyze uniqueness")
    print("A known theorem states that cr(G') = 1 for this type of graph.")
    print("The question is whether the 1-crossing drawing is unique.")
    print("Consider a planar drawing of G:")
    print("Draw vertices 0, 1, 2 as a triangle. Place vertex 3 inside this triangle, and vertex 4 outside.")
    print("Connect 3 and 4 to all of 0, 1, and 2.")
    print("\nTo add the edge e = (3, 4), the line must connect the inside vertex to the outside vertex.")
    print("This line MUST cross the boundary of the triangle (0, 1, 2).")
    print("There are three choices for which edge to cross:")
    print("1. Cross edge (0, 1)")
    print("2. Cross edge (1, 2)")
    print("3. Cross edge (0, 2)")
    print("\nSince there are multiple choices for which edge of G gets crossed by the new edge e,")
    print("the resulting 1-crossing drawing of G' is NOT unique.")
    print("-" * 30)
    
    print("Conclusion: The correct statement is that G' can be drawn with at most one crossing, but this drawing is not unique.")

demonstrate_maximal_planar_properties()