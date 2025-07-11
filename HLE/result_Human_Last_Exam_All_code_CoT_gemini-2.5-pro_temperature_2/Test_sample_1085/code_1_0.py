import collections

def solve_prime_path():
    """
    Calculates the number of distinct Prime Paths of 4 moves
    from (1, 1) to (5, 7) in PrimeGrid+1.
    """
    # Define the sequence of primes plus 1. A 4-step path won't go too far,
    # so we only need primes up to a certain point.
    P_sequence = [1, 2, 3, 5, 7, 11, 13]
    P_indices = {p: i for i, p in enumerate(P_sequence)}

    def get_neighbors(node):
        """Finds adjacent Prime Intersections for a given node."""
        x, y = node
        neighbors = []
        
        # Horizontal neighbors
        if x in P_indices:
            idx = P_indices[x]
            if idx > 0: # a predecessor exists
                neighbors.append((P_sequence[idx - 1], y))
            if idx < len(P_sequence) - 1: # a successor exists
                neighbors.append((P_sequence[idx + 1], y))

        # Vertical neighbors
        if y in P_indices:
            idx = P_indices[y]
            if idx > 0: # a predecessor exists
                neighbors.append((x, P_sequence[idx - 1]))
            if idx < len(P_sequence) - 1: # a successor exists
                neighbors.append((x, P_sequence[idx + 1]))
                
        return neighbors

    # Let N(node, k) be the number of distinct paths of length k from (1,1).
    # We will calculate this iteratively.
    start_node = (1, 1)
    end_node = (5, 7)
    num_moves = 4
    
    # path_trace stores N(node, k) for each step k.
    path_trace = {0: collections.defaultdict(int, {start_node: 1})}

    for step in range(1, num_moves + 1):
        paths_at_next_step = collections.defaultdict(int)
        # Get nodes from the previous step
        prev_paths = path_trace[step - 1]
        for node, count in prev_paths.items():
            for neighbor in get_neighbors(node):
                paths_at_next_step[neighbor] += count
        path_trace[step] = paths_at_next_step
        
    # Explain the logic and derive the final answer.
    print("To find the number of paths of length 4 to (5, 7), we must find the number of paths of length 3 to its neighbors.")
    
    end_node_neighbors = get_neighbors(end_node)
    print(f"The neighbors of the destination (5, 7) are {end_node_neighbors}.\n")
    
    # For each neighbor of the end_node, find the number of paths of length 3 to it.
    final_sum_components = []
    paths_to_neighbors_at_step_2 = path_trace[2]

    for neighbor_of_end_node in end_node_neighbors:
        # Paths to this neighbor at step 3 must come from its own neighbors at step 2.
        count_from_step2 = 0
        neighbors_of_neighbor = get_neighbors(neighbor_of_end_node)
        
        print(f"Calculating paths to {neighbor_of_end_node} (at step 3):")
        print(f"  - This requires paths to its neighbors {neighbors_of_neighbor} at step 2.")
        
        for n in neighbors_of_neighbor:
            count_from_step2 += paths_to_neighbors_at_step_2.get(n, 0)
        
        print(f"  - The number of paths from step 2 to these nodes is {count_from_step2}.")
        print(f"  - Therefore, N({neighbor_of_end_node}, 3) = {count_from_step2}.\n")
        final_sum_components.append(count_from_step2)

    final_result = sum(final_sum_components)
    
    # Format the final equation as requested.
    equation_str = " + ".join(map(str, final_sum_components))
    print("Final Calculation:")
    print(f"Total Paths = (Paths to (3, 7) at step 3) + (Paths to (5, 5) at step 3)")
    print(f"Total Paths = {equation_str} = {final_result}")

if __name__ == '__main__':
    solve_prime_path()