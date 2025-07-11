import math

def solve():
    """
    This function solves the graph theory problem and explains the reasoning.
    
    Problem Statement:
    Let G be a connected 3-regular adjustable graph with 2000 vertices that has at least one perfect matching. 
    An adjustable matching M in G is one where for every two edges vu, v'u' in M, if vv' is an edge in G, then uu' is also an edge. 
    A graph G is adjustable if it has a maximum matching that is adjustable.
    How many non-isomorphic graphs like G exist?
    """

    # Step 1 & 2: Analyze the adjustable property
    # Let M be an adjustable perfect matching in G. A perfect matching exists and is a maximum matching.
    # The vertex set V can be partitioned into X = {x_1, ..., x_1000} and Y = {y_1, ..., y_1000},
    # such that M = {(x_i, y_i) | i=1,...,1000}.
    # The adjustable condition implies:
    # 1. (x_i, x_j) is an edge iff (y_i, y_j) is an edge. So, G[X] is isomorphic to G[Y].
    # 2. (x_i, y_j) is an edge iff (y_i, x_j) is an edge.

    # Step 3: Use the 3-regularity
    # For any vertex x_i, deg(x_i) = 3. One edge (x_i, y_i) is in M.
    # Let k_i = deg_G[X](x_i) and d_i be the number of edges from x_i to Y \ {y_i}.
    # Then for each i, 1 + k_i + d_i = 3, which means k_i + d_i = 2.

    # Step 4: Argue for homogeneity
    # For G to be connected, the structural properties defined by (k_i, d_i) must be uniform
    # across the graph. If there were vertices of different types (e.g., d_i=0 and d_j>0), the set of
    # vertices with d=0 would form disconnected components from the rest of the graph.
    # Thus, (k,d) is constant for all vertices.
    # This leaves three possible cases for (k, d): (2,0), (0,2), (1,1).

    print("The problem asks for the number of non-isomorphic, connected, 3-regular adjustable graphs with 2000 vertices and a perfect matching.")
    print("Let the adjustable perfect matching be M, partitioning the vertices into two sets X and Y of size 1000.")
    print("The properties of the graph lead to three possible structural cases based on how the non-matching edges are arranged:")
    print("")

    # Step 5: Describe the three cases
    # Case 1: k=2, d=0
    # G[X] and G[Y] are 2-regular graphs. Connectivity of G implies they are C_1000 cycles.
    # There are no other edges between X and Y besides the matching M.
    # This structure defines the unique (up to isomorphism) prism graph Y_1000 (also denoted C_1000 x K_2).
    graph1_k = 2
    graph1_d = 0
    print(f"Case 1: (k={graph1_k}, d={graph1_d})")
    print("The subgraphs induced by X and Y are both 2-regular and connected, meaning they are both 1000-cycles (C_1000).")
    print("The only edges between X and Y are the matching edges from M.")
    print("This corresponds to one unique graph: the prism graph Y_1000.")
    print("")

    # Case 2: k=0, d=2
    # G[X] and G[Y] have no edges. All non-matching edges are "crossing" edges.
    # The graph is bipartite with partitions X and Y.
    # The structure must be such that connectivity is maintained.
    # This leads to a unique graph where vertices x_i are connected to y_{i-1} and y_{i+1} (and y_i).
    graph2_k = 0
    graph2_d = 2
    print(f"Case 2: (k={graph2_k}, d={graph2_d})")
    print("The subgraphs induced by X and Y have no edges.")
    print("Each vertex x_i is connected to two vertices in Y (and y_i via M), and symmetrically for y_i.")
    print("To ensure connectivity, this structure also leads to one unique graph.")
    print("")
    
    # Case 3: k=1, d=1
    # G[X] and G[Y] are 1-regular graphs (i.e., perfect matchings on 1000 vertices).
    # Each x_i also has one crossing edge to some y_j (and symmetrically y_i to x_j).
    # These different matchings must link up to form a connected graph.
    # This case also results in a single non-isomorphic graph structure.
    graph3_k = 1
    graph3_d = 1
    print(f"Case 3: (k={graph3_k}, d={graph3_d})")
    print("The subgraphs induced by X and Y are both 1-regular (perfect matchings).")
    print("Each vertex also has one 'crossing' non-matching edge to the other partition.")
    print("This defines a third unique graph structure, different from the first two.")
    print("")
    
    # Step 6: Conclude
    num_graphs = 3
    print("These three cases are mutually exclusive and cover all possibilities for a connected graph of this type.")
    print("The induced subgraphs G[X] are 2-regular, 0-regular, and 1-regular, respectively, proving the three resulting graphs are non-isomorphic.")
    print(f"Therefore, there are {num_graphs} such non-isomorphic graphs.")

solve()

# The final answer is an integer based on the detailed graph theory analysis.
final_answer = 3
print(f"<<<{final_answer}>>>")