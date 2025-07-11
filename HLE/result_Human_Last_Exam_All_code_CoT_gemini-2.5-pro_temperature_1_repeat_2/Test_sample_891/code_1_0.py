from collections import deque

def solve_oil_division():
    """
    Solves the oil division puzzle using Breadth-First Search (BFS)
    to find the shortest sequence of pouring actions.
    """
    # Capacities: (X is inf, A=8, B=17, C=21)
    # Using tuple for state: (oil_in_x, oil_in_a, oil_in_b, oil_in_c)
    capacities = (float('inf'), 8, 17, 21)
    container_names = ['X', 'A', 'B', 'C']
    
    initial_state = (39, 0, 0, 0)
    target_state = (13, 0, 13, 13)

    # The queue will store tuples of (state, path_to_state)
    queue = deque([(initial_state, [])])
    
    # The visited set stores states we've already processed to avoid cycles
    visited = {initial_state}

    while queue:
        current_state, path = queue.popleft()

        if current_state == target_state:
            # Solution found, print the results
            print("Shortest sequence of operations found:")
            
            # Trace and print each step
            s = initial_state
            print(f"Start: {s}")
            for i, action in enumerate(path):
                source_idx, dest_idx = action
                source_amount = s[source_idx]
                dest_space = capacities[dest_idx] - s[dest_idx]
                pour = min(source_amount, dest_space)
                
                new_s_list = list(s)
                new_s_list[source_idx] -= pour
                new_s_list[dest_idx] += pour
                s = tuple(int(v) for v in new_s_list)
                
                print(f"Step {i+1}: P({container_names[source_idx]}, {container_names[dest_idx]}) -> {s}")

            # Print the final answer as requested
            final_op_source, final_op_dest = path[-1]
            print("\n---")
            print(f"The final state is X={s[0]}, A={s[1]}, B={s[2]}, C={s[3]}.")
            print(f"The three equal portions are in containers X, B, and C.")
            print(f"Final equation: {s[0]} = {s[2]} = {s[3]}")
            print(f"The last operation is P({container_names[final_op_source]}, {container_names[final_op_dest]}).")
            
            # Determine the answer choice letter
            choices = { "A": "P(A, B)", "B": "P(A, C)", "C": "P(B, C)", "D": "P(C, B)",
                        "E": "P(C, A)", "F": "P(X, A)", "G": "P(X, B)", "H": "P(B, A)" }
            last_op_str = f"P({container_names[final_op_source]}, {container_names[final_op_dest]})"
            answer_letter = "Unknown"
            for letter, op_str in choices.items():
                if op_str == last_op_str:
                    answer_letter = letter
                    break
            print(f"This corresponds to answer choice: {answer_letter}")
            return f"<<<{answer_letter}>>>"

        # Generate all possible next states
        # indices: 0=X, 1=A, 2=B, 3=C
        for source_idx in range(4):
            # Destination cannot be X (index 0)
            for dest_idx in range(1, 4):
                if source_idx == dest_idx:
                    continue
                
                # Check if a pour is possible
                source_amount = current_state[source_idx]
                dest_amount = current_state[dest_idx]
                dest_capacity = capacities[dest_idx]
                
                if source_amount > 0 and dest_amount < dest_capacity:
                    pour_amount = min(source_amount, dest_capacity - dest_amount)
                    
                    new_state_list = list(current_state)
                    new_state_list[source_idx] -= pour_amount
                    new_state_list[dest_idx] += pour_amount
                    new_state = tuple(new_state_list)

                    if new_state not in visited:
                        visited.add(new_state)
                        new_path = path + [(source_idx, dest_idx)]
                        queue.append((new_state, new_path))
                        
    return "No solution found."

# Execute the solver
solve_oil_division()