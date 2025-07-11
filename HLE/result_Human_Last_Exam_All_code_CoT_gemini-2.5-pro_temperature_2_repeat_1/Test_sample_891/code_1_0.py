from collections import deque

def solve_pouring_problem():
    """
    Solves the oil division puzzle using Breadth-First Search (BFS)
    to find the shortest sequence of pouring actions.
    """
    # Capacities are ignored for X since we can't pour back into it.
    capacities = {'X': 39, 'A': 8, 'B': 17, 'C': 21}
    container_names = ['X', 'A', 'B', 'C']
    
    # Initial state (X, A, B, C)
    start_state = (39, 0, 0, 0)
    
    # Target state: 13L in X, B, and C, leaving 0L in A.
    target_state = (13, 0, 13, 13)

    # Queue for BFS: stores tuples of (state_tuple, path_list)
    queue = deque([(start_state, [])])
    
    # Set to keep track of visited states to avoid cycles.
    visited = {start_state}
    
    while queue:
        current_state, path = queue.popleft()

        # Check if the target state is reached
        if current_state == target_state:
            print("Shortest sequence found:\n")
            
            # Print the sequence of operations with state changes
            state = list(start_state)
            print(f"Start: X={state[0]}, A={state[1]}, B={state[2]}, C={state[3]}")
            
            for i, move in enumerate(path):
                source_name, dest_name = move.split('->')
                s_idx, d_idx = container_names.index(source_name), container_names.index(dest_name)
                
                source_val = state[s_idx]
                dest_val = state[d_idx]
                dest_cap = capacities[dest_name]

                amount_to_pour = min(source_val, dest_cap - dest_val)
                
                state[s_idx] -= amount_to_pour
                state[d_idx] += amount_to_pour
                
                print(f"{i + 1}. P({source_name}, {dest_name}): X={state[0]}, A={state[1]}, B={state[2]}, C={state[3]}")

            final_state = state
            print("\nFinal State Achieved: 13L in X, 13L in B, and 13L in C.")
            print("Final Equation:", f"{final_state[0]} + {final_state[1]} + {final_state[2]} + {final_state[3]} = 39")
            print("In three portions:", f"{final_state[0]} + {final_state[2]} + {final_state[3]} = 39")
            
            # Determine the choice letter for the last operation
            last_op_str = path[-1].replace('->', ', ')
            last_op_formatted = f"P({last_op_str})"
            
            # Map the final action to the provided choices
            choice_map = {
                "P(A, B)": "A", "P(A, C)": "B", "P(B, C)": "C",
                "P(C, B)": "D", "P(C, A)": "E", "P(X, A)": "F",
                "P(X, B)": "G", "P(B, A)": "H"
            }
            final_choice = choice_map.get(last_op_formatted, "Unknown")
            
            print(f"\nThe last operation is {last_op_formatted}.")
            print(f"<<<{final_choice}>>>")
            return

        # Generate next possible states
        current_state_list = list(current_state)
        for i in range(len(container_names)):  # Source container index
            for j in range(len(container_names)): # Destination container index
                if i == j:
                    continue

                source_name = container_names[i]
                dest_name = container_names[j]

                # Restriction: cannot pour back to container X
                if dest_name == 'X':
                    continue

                source_val = current_state_list[i]
                dest_val = current_state_list[j]
                dest_cap = capacities[dest_name]
                
                # If source is empty or destination is full, can't pour
                if source_val == 0 or dest_val == dest_cap:
                    continue

                # This single calculation handles both "empty source" and "fill dest" rules
                amount_to_pour = min(source_val, dest_cap - dest_val)
                
                next_state_list = list(current_state)
                next_state_list[i] -= amount_to_pour
                next_state_list[j] += amount_to_pour
                
                next_state = tuple(next_state_list)

                if next_state not in visited:
                    visited.add(next_state)
                    new_path = path + [f"{source_name}->{dest_name}"]
                    queue.append((next_state, new_path))
    
    print("No solution found.")

if __name__ == '__main__':
    solve_pouring_problem()