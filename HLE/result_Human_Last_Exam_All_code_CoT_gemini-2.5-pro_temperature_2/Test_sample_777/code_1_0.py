import networkx as nx
from itertools import combinations

def find_k_disjoint_long_cycles(graph, k):
    """
    Finds k vertex-disjoint cycles of length at least k using backtracking.
    This is a demonstration and not an efficient FPT algorithm.
    """
    
    # Pre-check: not enough vertices for k cycles of length k
    if len(graph.nodes) < k * k:
        return None
        
    # The main recursive backtracking function
    def backtrack_find(current_graph, k_rem):
        if k_rem == 0:
            return []  # Base case: successfully found all k cycles

        # In a real backtracking algorithm, we should iterate through all possible cycles.
        # nx.simple_cycles can be very slow. We'll find cycles from each node as a starting point.
        # This is not guaranteed to find ALL simple cycles but serves for demonstration.
        
        nodes_to_try = list(current_graph.nodes())
        
        for u in nodes_to_try:
            # We need to find simple cycles starting from u of length at least k
            # Using cycle_basis can be an alternative, but doesn't guarantee finding all simple cycles
            # For this example, let's use a robust (but potentially slow) cycle finding method
            # on subgraphs to explore different cycle choices.
            
            # The cycles must be found in the current_graph.
            for cycle in nx.simple_cycles(current_graph):
                if len(cycle) < k:
                    continue

                # Found a candidate cycle. Let's try to proceed with this choice.
                
                sub_solution = [cycle]
                
                # Create a new graph without the vertices of the found cycle
                graph_after_removal = current_graph.copy()
                graph_after_removal.remove_nodes_from(cycle)
                
                # Recurse to find the remaining k-1 cycles
                remaining_solution = backtrack_find(graph_after_removal, k_rem - 1)
                
                if remaining_solution is not None:
                    # Success! Found a valid solution down this path.
                    return sub_solution + remaining_solution

        # If we have tried all possibilities from this state and none led to a solution
        return None

    solution = backtrack_find(graph, k)
    return solution


if __name__ == '__main__':
    # --- Example Usage ---
    # Create a graph with a known solution
    G = nx.Graph()
    # Add two disjoint cycles, C5 and C4. Let k=2. Lengths must be >= 2.
    G.add_edges_from([(1, 2), (2, 3), (3, 4), (4, 5), (5, 1)])  # Cycle 1 (length 5)
    G.add_edges_from([(6, 7), (7, 8), (8, 9), (9, 6)])      # Cycle 2 (length 4)
    # Add one more vertex to make it more interesting
    G.add_edge(1, 10)

    k_param = 2
    
    print(f"Searching for k={k_param} vertex-disjoint cycles of length at least {k_param}...")

    disjoint_cycles = find_k_disjoint_long_cycles(G, k_param)

    if disjoint_cycles:
        print("\nFound a solution:")
        # The prompt mentions "output each number in the final equation".
        # This might refer to printing the parameters of the solution, which we do here.
        # For k=2, we found 2 cycles. Their lengths must be >= 2.
        total_vertices = 0
        for i, cycle in enumerate(disjoint_cycles):
            length = len(cycle)
            total_vertices += length
            # Here we output the numbers related to each part of the solution.
            print(f"  Cycle {i+1} (length = {length}): {cycle}")
        
        # We can also output the summary equation, as mentioned.
        # This confirms that the total number of vertices respects the k*k bound.
        print(f"\nFinal check: sum of cycle lengths = {total_vertices}, which is >= k*k = {k_param*k_param}.")
    else:
        print("\nNo solution found.")
        print("The graph does not contain k disjoint cycles of the required length.")