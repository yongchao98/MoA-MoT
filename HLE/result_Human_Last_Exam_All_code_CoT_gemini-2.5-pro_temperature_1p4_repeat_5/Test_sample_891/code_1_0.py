import collections

def solve_oil_division():
    """
    Solves the oil division puzzle using Breadth-First Search (BFS)
    to find the shortest sequence of pouring actions.
    """
    # Define container names and their capacities.
    # The initial amount in X is treated as its capacity for state representation.
    capacities = {'X': 39, 'A': 8, 'B': 17, 'C': 21}
    container_indices = {'X': 0, 'A': 1, 'B': 2, 'C': 3}
    index_to_name = {v: k for k, v in container_indices.items()}
    
    # Define the initial and target states as tuples (X, A, B, C)
    initial_state = (39, 0, 0, 0)
    target_state = (13, 0, 13, 13)

    # The queue for BFS stores tuples of (state, path)
    queue = collections.deque([(initial_state, [])])
    
    # A set to keep track of visited states
    visited = {initial_state}

    # Begin the BFS search
    while queue:
        current_state_tuple, path = queue.popleft()

        # Check if the current state is the target
        if current_state_tuple == target_state:
            print("Solution Found!")
            print(f"The shortest sequence requires {len(path)} steps.")
            print("-" * 30)
            
            # Print the full sequence of operations
            temp_state = initial_state
            print(f"Start State: (X={temp_state[0]}, A={temp_state[1]}, B={temp_state[2]}, C={temp_state[3]})")
            for i, move in enumerate(path):
                # To display the state after each move, we recalculate it here
                source_name, dest_name = move.replace('P(', '').replace(')', '').split(', ')
                s_idx, d_idx = container_indices[source_name], container_indices[dest_name]
                
                state_list = list(temp_state)
                amount = min(state_list[s_idx], capacities[dest_name] - state_list[d_idx])
                state_list[s_idx] -= amount
                state_list[d_idx] += amount
                temp_state = tuple(state_list)
                
                print(f"Step {i+1}: {move} -> State: (X={temp_state[0]}, A={temp_state[1]}, B={temp_state[2]}, C={temp_state[3]})")

            print("-" * 30)
            print(f"Final action to achieve the goal: {path[-1]}")
            return path[-1]

        # Generate all valid next states by trying every possible pour
        # 'i' is the source container index, 'j' is the destination container index
        current_state_list = list(current_state_tuple)
        for i in range(4):  # Source: X, A, B, C
            for j in range(1, 4):  # Destination: A, B, C (cannot pour to X)
                if i == j:
                    continue

                source_name = index_to_name[i]
                dest_name = index_to_name[j]

                # If source is empty or destination is full, skip
                if current_state_list[i] == 0 or current_state_list[j] == capacities[dest_name]:
                    continue
                
                next_state_list = list(current_state_list)
                
                # Calculate amount to pour
                amount_to_pour = min(next_state_list[i], capacities[dest_name] - next_state_list[j])
                
                # Perform the pour
                next_state_list[i] -= amount_to_pour
                next_state_list[j] += amount_to_pour
                
                next_state_tuple = tuple(next_state_list)
                
                if next_state_tuple not in visited:
                    visited.add(next_state_tuple)
                    new_path = path + [f"P({source_name}, {dest_name})"]
                    queue.append((next_state_tuple, new_path))
                    
    print("No solution was found.")
    return None

# Execute the solver function
solve_oil_division()