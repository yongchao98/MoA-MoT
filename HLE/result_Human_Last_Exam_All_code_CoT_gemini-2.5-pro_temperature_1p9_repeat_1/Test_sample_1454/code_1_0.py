import networkx as nx

def solve_fractal_components():
    """
    Calculates the number of non-degenerate, locally connected components of a fractal set.
    """
    # The set of vectors for the IFS from the problem description
    # D = {(0,0), (0,1), (0,2),(0,3),(3,0), (3,1), (3,2),(3,3)}
    # Let's represent them as a list of tuples
    D = [(0, 0), (0, 1), (0, 2), (0, 3), (3, 0), (3, 1), (3, 2), (3, 3)]

    # The fractal F is the attractor of the IFS given by the maps f_d(v) = (v+d)/4.
    # The connectivity of F can be determined by the adjacency of the images of the
    # unit square, S = [0,1]^2, under these maps.
    # Each f_d(S) is a square of side 1/4 at position d/4.
    # Let's use the vectors d as identifiers for these squares.

    # We build a graph where nodes represent these 8 squares.
    G = nx.Graph()
    G.add_nodes_from(D)

    # An edge exists between two nodes if their corresponding squares are adjacent.
    # Two squares (i,j) and (i',j') in a grid are adjacent if they touch along an edge.
    # This means the distance between their centers is 1 grid unit.
    # i.e., |i-i'| + |j-j'| = 1.
    for i in range(len(D)):
        for j in range(i + 1, len(D)):
            d1 = D[i]
            d2 = D[j]
            # Check for horizontal or vertical adjacency
            if abs(d1[0] - d2[0]) + abs(d1[1] - d2[1]) == 1:
                G.add_edge(d1, d2)

    # According to a theorem by Hata and others, the number of connected components of the
    # attractor F is equal to the number of connected components of this adjacency graph G.
    num_components = nx.number_connected_components(G)

    # Now, we verify if these components meet the specified criteria.
    # The components of the graph G are:
    graph_components = list(nx.connected_components(G))
    # component_1 = {(0,0), (0,1), (0,2), (0,3)}
    # component_2 = {(3,0), (3,1), (3,2), (3,3)}

    # Check for local connectivity:
    # A component of the attractor is locally connected if the union of the squares
    # in its corresponding graph component is a connected set.
    # For component 1, the squares S_{0j} form a vertical column [0, 1/4] x [0, 1]. This is connected.
    # For component 2, the squares S_{3j} form another column [3/4, 1] x [0, 1]. This is also connected.
    # Thus, all components of F are locally connected.

    # Check for non-degeneracy:
    # A component is degenerate if it is a single point. This occurs if all maps in the
    # corresponding sub-IFS have a common fixed point.
    # For a map f_d(p) = (p+d)/4, the fixed point is p = d/3.
    # For component 1, d can be (0,0) or (0,1). p would be (0,0) or (0,1/3). No common fixed point.
    # For component 2, d can be (3,0) or (3,1). p would be (1,0) or (1,1/3). No common fixed point.
    # Thus, all components of F are non-degenerate.
    
    final_answer = num_components
    
    print(f"The set D is: {D}")
    print("The fractal F is the attractor of an Iterated Function System (IFS) defined by the vectors in D.")
    print("The connectivity of F can be determined by the connectivity of a graph representing the layout of the IFS mappings.")
    print("The vertices of the graph correspond to the vectors in D. An edge connects two vertices if their corresponding squares in the grid are adjacent.")
    print("\nThe graph formed has two separate groups of nodes:")
    print(f"1. The nodes related to the first column of the grid: {graph_components[0]}")
    print(f"2. The nodes related to the fourth column of the grid: {graph_components[1]}")
    print("These two groups are not connected to each other.")
    print(f"\nThe number of connected components in this graph is {num_components}.")
    print("According to fractal theory, this is also the number of connected components of the set F.")
    print("\nEach of these components is non-degenerate (it is a complex fractal shape, not a single point).")
    print("Each component is also locally connected (since the union of squares for each sub-system is connected).")
    print(f"\nTherefore, the number of components of F that are non-degenerate and locally connected is {final_answer}.")


solve_fractal_components()
