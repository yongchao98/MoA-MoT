import collections

def solve_pouring_puzzle():
    """
    Solves the oil pouring puzzle using Breadth-First Search to find the shortest sequence.
    """
    # Capacities of the containers A, B, and C.
    capacities = {'X': float('inf'), 'A': 8, 'B': 17, 'C': 21}
    
    # Define the starting state and the target state.
    # State is a tuple: (oil_in_X, oil_in_A, oil_in_B, oil_in_C)
    initial_state = (39, 0, 0, 0)
    target_state = (13, 0, 13, 13)

    # The queue for BFS will store tuples of (state, path_to_state).
    queue = collections.deque([(initial_state, [])])
    
    # A set to keep track of visited states to avoid cycles and redundant work.
    visited = {initial_state}

    while queue:
        current_state, path = queue.popleft()

        # If the target is reached, we found the shortest path.
        if current_state == target_state:
            # The puzzle is solved, print the final state and the last operation.
            final_x, final_a, final_b, final_c = current_state
            print(f"Final state reached: X={final_x}, A={final_a}, B={final_b}, C={final_c}")
            last_operation = path[-1]
            print(f"The last operation is: {last_operation}")
            return path

        # Generate all possible next states from the current state.
        state_map = {'X': current_state[0], 'A': current_state[1], 'B': current_state[2], 'C': current_state[3]}
        
        # All possible pouring combinations (source, destination).
        all_pours = [('X', 'A'), ('X', 'B'), ('X', 'C'),
                     ('A', 'B'), ('A', 'C'),
                     ('B', 'A'), ('B', 'C'),
                     ('C', 'A'), ('C', 'B')]

        for s_name, d_name in all_pours:
            s_vol = state_map[s_name]
            d_vol = state_map[d_name]
            d_cap = capacities[d_name]

            # Cannot pour from an empty container or into a full one.
            if s_vol == 0 or d_vol == d_cap:
                continue

            # Case 1: Fill the destination container.
            # This action is possible if the source has enough liquid to fill the destination.
            amount_to_fill = d_cap - d_vol
            if s_vol >= amount_to_fill:
                next_state_map = state_map.copy()
                next_state_map[s_name] -= amount_to_fill
                next_state_map[d_name] += amount_to_fill
                next_state_tuple = (next_state_map['X'], next_state_map['A'], next_state_map['B'], next_state_map['C'])
                
                if next_state_tuple not in visited:
                    visited.add(next_state_tuple)
                    new_path = path + [f"P({s_name}, {d_name})"]
                    queue.append((next_state_tuple, new_path))
            
            # Case 2: Empty the source container.
            # This action is possible if the destination has enough space for all the liquid from the source.
            amount_to_empty = s_vol
            if d_cap - d_vol >= amount_to_empty:
                next_state_map = state_map.copy()
                next_state_map[s_name] -= amount_to_empty
                next_state_map[d_name] += amount_to_empty
                next_state_tuple = (next_state_map['X'], next_state_map['A'], next_state_map['B'], next_state_map['C'])
                
                # A single pour can satisfy both conditions (e.g., empty A fills C exactly).
                # We check `visited` again to avoid adding the same state twice in such cases.
                if next_state_tuple not in visited:
                    visited.add(next_state_tuple)
                    new_path = path + [f"P({s_name}, {d_name})"]
                    queue.append((next_state_tuple, new_path))
                    
    return None # Should not be reached if a solution exists

if __name__ == '__main__':
    solution_path = solve_pouring_puzzle()
    if not solution_path:
        print("No solution was found.")
