import networkx as nx
import numpy as np

def solve_graph_transformation_problem():
    """
    Finds the number of non-isomorphic, connected graphs G with at most 7 vertices
    that satisfy T(G) = G, where T is the specified graph transformation.
    The script prints the counts for each number of vertices and the total sum.
    """
    # Use a dictionary to store counts for each number of vertices n=1 to 7.
    counts_by_n = {n: 0 for n in range(1, 8)}

    # networkx.graph_atlas_g() returns a list of all non-isomorphic graphs
    # with up to 7 vertices.
    graph_list = nx.graph_atlas_g()

    for G in graph_list:
        # Per the problem, we only consider connected graphs.
        if not nx.is_connected(G):
            continue
            
        n = G.number_of_nodes()

        # The atlas contains graphs from n=0 to n=7. We only need 1 to 7.
        if n == 0:
            continue

        # Get the adjacency matrix A for the graph G.
        # Sorting nodes ensures a consistent matrix representation.
        nodes = sorted(G.nodes())
        A = nx.to_numpy_array(G, nodelist=nodes)

        # The number of length-2 paths between any two vertices is given by A^2.
        A_squared = np.dot(A, A)

        # Now, construct the adjacency matrix for the transformed graph T(G).
        # An edge (x,y) exists in T(G) iff there are 1 or 2 length-2 paths in G.
        # This translates to A_transformed[i,j] = 1 iff 1 <= A_squared[i,j] <= 2.
        A_transformed = ((A_squared >= 1) & (A_squared <= 2)).astype(np.int64)

        # The transformation T does not allow self-loops, so the diagonal must be zero.
        np.fill_diagonal(A_transformed, 0)

        # We are looking for graphs where T(G) = G. This means their adjacency
        # matrices must be identical.
        if np.array_equal(A, A_transformed):
            if n in counts_by_n:
                counts_by_n[n] += 1

    # Format the output to show the count for each vertex number and the total sum,
    # satisfying the instruction "output each number in the final equation".
    equation_parts = []
    for n in sorted(counts_by_n.keys()):
        equation_parts.append(str(counts_by_n[n]))
    
    total_count = sum(counts_by_n.values())
    equation_str = " + ".join(equation_parts)
    
    print("Found the following number of graphs for n = 1, 2, 3, 4, 5, 6, and 7 respectively:")
    print(f"Final Count = {equation_str} = {total_count}")

# Execute the function to find the answer.
solve_graph_transformation_problem()