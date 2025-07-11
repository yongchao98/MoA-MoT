from collections import deque

def solve_oil_division_problem():
    """
    Solves the oil division puzzle using Breadth-First Search (BFS)
    to find the shortest sequence of pouring actions.
    """
    # Define container capacities and names
    capacities = {'X': float('inf'), 'A': 8, 'B': 17, 'C': 21}
    names = ['X', 'A', 'B', 'C']
    name_to_idx = {name: i for i, name in enumerate(names)}

    # Define initial and goal states
    initial_state = (39, 0, 0, 0)
    goal_state = (0, 13, 13, 13)

    # All possible pouring actions (source_name, destination_name)
    # Pouring back to 'X' is not allowed.
    possible_moves = [
        ('X', 'A'), ('X', 'B'), ('X', 'C'),
        ('A', 'B'), ('A', 'C'),
        ('B', 'A'), ('B', 'C'),
        ('C', 'A'), ('C', 'B'),
    ]

    # Queue for BFS: stores tuples of (state_tuple, path_list)
    queue = deque([(initial_state, [])])
    # Set to keep track of visited states to avoid cycles
    visited = {initial_state}

    # Helper function to perform a pour and return the new state
    def perform_pour(state, source_name, dest_name):
        state_list = list(state)
        s_idx, d_idx = name_to_idx[source_name], name_to_idx[dest_name]
        
        source_val = state_list[s_idx]
        dest_val = state_list[d_idx]
        dest_capacity = capacities[dest_name]

        amount_to_pour = min(source_val, dest_capacity - dest_val)

        if amount_to_pour == 0:
            return None # No change in state

        state_list[s_idx] -= amount_to_pour
        state_list[d_idx] += amount_to_pour
        return tuple(state_list)

    solution_path = None
    while queue:
        current_state, path = queue.popleft()

        if current_state == goal_state:
            solution_path = path
            break

        for source_name, dest_name in possible_moves:
            new_state = perform_pour(current_state, source_name, dest_name)
            
            if new_state and new_state not in visited:
                visited.add(new_state)
                new_path = path + [f"P({source_name}, {dest_name})"]
                queue.append((new_state, new_path))
    
    # --- Output the results ---
    if solution_path:
        print("Solution Found! Here is the shortest sequence of operations:\n")
        
        # Print the trace of actions and resulting states
        state = initial_state
        print(f"Start: X={state[0]:>2}, A={state[1]:>2}, B={state[2]:>2}, C={state[3]:>2}")
        
        for i, action_str in enumerate(solution_path):
            # Parse action string to get source and destination
            source, dest = action_str[2], action_str[5]
            state = perform_pour(state, source, dest)
            print(f"{i+1}. {action_str}: X={state[0]:>2}, A={state[1]:>2}, B={state[2]:>2}, C={state[3]:>2}")

        # The question is about the last operation
        last_operation = solution_path[-1]
        print(f"\nThe last operation to achieve the goal is: {last_operation}")

        # Map the answer to the provided choices
        choices = {
            'A': 'P(A, B)', 'B': 'P(A, C)', 'C': 'P(B, C)', 'D': 'P(C, B)',
            'E': 'P(C, A)', 'F': 'P(X, A)', 'G': 'P(X, B)', 'H': 'P(B, A)'
        }
        
        final_answer_key = None
        for key, value in choices.items():
            if value == last_operation:
                final_answer_key = key
                break
        
        if final_answer_key:
            print(f"<<<{final_answer_key}>>>")
        else:
            print("Error: The found operation is not in the answer choices.")
    else:
        print("No solution was found.")

# Execute the solver
solve_oil_division_problem()