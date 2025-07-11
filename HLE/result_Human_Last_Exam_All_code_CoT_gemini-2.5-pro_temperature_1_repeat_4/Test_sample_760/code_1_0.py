import networkx as nx
import numpy as np

def check_graph_is_fixed_point(G):
    """
    Checks if a graph G satisfies the transformation T(G) = G.

    The transformation T adds an edge between distinct vertices x and y
    if and only if there are exactly one or two length-2 paths between them.

    Args:
        G: A networkx Graph object.

    Returns:
        True if T(G) = G, False otherwise.
    """
    n = G.number_of_nodes()
    if n <= 1:
        # The single-vertex graph is connected and T(G)=G.
        return True

    # Get the adjacency matrix A of the graph G.
    # Sorting nodes ensures a consistent matrix representation.
    nodelist = sorted(G.nodes())
    A = nx.to_numpy_array(G, nodelist=nodelist)

    # A^2 gives the number of length-2 paths between vertices.
    A_squared = np.dot(A, A)

    # Construct the adjacency matrix A_prime for the transformed graph T(G).
    A_prime = np.zeros_like(A)
    for i in range(n):
        for j in range(i + 1, n):
            # The condition for an edge in T(G) is 1 or 2 paths of length 2.
            if 1 <= A_squared[i, j] <= 2:
                A_prime[i, j] = 1
                A_prime[j, i] = 1

    # Check if the original graph G is identical to the transformed graph T(G).
    return np.array_equal(A, A_prime)

def solve_problem():
    """
    Finds the number of non-isomorphic, connected graphs with at most 7 vertices
    that are fixed points of the transformation T.
    """
    counts_per_n = []
    
    # Iterate through number of vertices from 1 to 7.
    for n in range(1, 8):
        count_for_n = 0
        # Generate all non-isomorphic graphs with n vertices.
        graph_iterator = nx.nonisomorphic_graphs(n)
        
        for G in graph_iterator:
            # We are only interested in connected graphs.
            if nx.is_connected(G):
                # Check if the graph is a fixed point of the transformation T.
                if check_graph_is_fixed_point(G):
                    count_for_n += 1
        
        counts_per_n.append(count_for_n)

    # Calculate the total count.
    total_count = sum(counts_per_n)

    # Format the output string as an equation.
    equation_parts = []
    for i, count in enumerate(counts_per_n):
        n = i + 1
        equation_parts.append(f"{count}")
    
    equation_str = " + ".join(equation_parts)
    
    print(f"The number of graphs for n=1 to 7 are: {equation_str}")
    print(f"The final equation is: {equation_str} = {total_count}")


if __name__ == "__main__":
    solve_problem()
