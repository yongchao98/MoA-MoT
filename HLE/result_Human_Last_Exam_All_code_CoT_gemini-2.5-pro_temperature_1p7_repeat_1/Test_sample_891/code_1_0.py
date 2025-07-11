import collections

def solve_oil_division():
    """
    Solves the oil division puzzle using Breadth-First Search (BFS)
    to find the shortest sequence of pouring actions.
    """
    
    # Define capacities of containers A, B, C. X's capacity is effectively infinite.
    # The tuple is indexed as (X, A, B, C)
    capacities = (float('inf'), 8, 17, 21)
    
    # Define the starting and goal states of the oil distribution.
    start_state = (39, 0, 0, 0)
    goal_state = (13, 0, 13, 13)
    
    # Initialize a queue for BFS. Each element is a tuple: (state, path_to_state)
    queue = collections.deque([(start_state, [])])
    
    # A set to store visited states to avoid cycles and redundant work.
    visited = {start_state}
    
    # Helper dictionary to map container indices to names for printing.
    container_names = {0: "X", 1: "A", 2: "B", 3: "C"}
    
    # Start the BFS loop.
    while queue:
        current_state, path = queue.popleft()
        
        # Check if the current state is the goal state.
        if current_state == goal_state:
            print("Shortest sequence to achieve the goal found:")
            
            # Print the step-by-step process for verification.
            state_tracker = list(start_state)
            print(f"Start: X={state_tracker[0]:2d}, A={state_tracker[1]:2d}, B={state_tracker[2]:2d}, C={state_tracker[3]:2d}")
            
            for i, op in enumerate(path):
                # Deconstruct operation string like "P(X,A)" to find source and destination
                source_name = op[2]
                dest_name = op[5]
                
                source_idx = [k for k, v in container_names.items() if v == source_name][0]
                dest_idx = [k for k, v in container_names.items() if v == dest_name][0]

                # Calculate the amount to pour based on the rules.
                amount_to_pour = min(state_tracker[source_idx], capacities[dest_idx] - state_tracker[dest_idx])
                
                # Update the state tracker.
                state_tracker[source_idx] -= amount_to_pour
                state_tracker[dest_idx] += amount_to_pour
                
                print(f"Step {i+1}: {op} -> X={state_tracker[0]:2d}, A={state_tracker[1]:2d}, B={state_tracker[2]:2d}, C={state_tracker[3]:2d}")

            # Print the final numbers in the goal state.
            print("\nFinal State Achieved:")
            print(f"Container X: {current_state[0]} liters")
            print(f"Container A: {current_state[1]} liters")
            print(f"Container B: {current_state[2]} liters")
            print(f"Container C: {current_state[3]} liters")
            print(f"Total: {current_state[0]} + {current_state[1]} + {current_state[2]} + {current_state[3]} = {sum(current_state)}")
            
            # Output the last operation, which is the answer to the question.
            last_operation = path[-1]
            print(f"\nThe last operation is: {last_operation}")

            # Map the operation string to the answer choices.
            if last_operation == "P(A, B)": answer = "A"
            elif last_operation == "P(A, C)": answer = "B"
            elif last_operation == "P(B, C)": answer = "C"
            elif last_operation == "P(C, B)": answer = "D"
            elif last_operation == "P(C, A)": answer = "E"
            elif last_operation == "P(X, A)": answer = "F"
            elif last_operation == "P(X, B)": answer = "G"
            elif last_operation == "P(B, A)": answer = "H"
            else: answer = "Unknown"
            
            print(f'<<<{answer}>>>')
            return

        # If not the goal, generate all possible next states from the current state.
        for i in range(4):  # Source container index (0:X, 1:A, 2:B, 3:C)
            for j in range(1, 4):  # Destination container index (1:A, 2:B, 3:C)
                
                # Can't pour from a container to itself.
                if i == j:
                    continue
                
                # Create a mutable list from the current state tuple.
                next_state_list = list(current_state)
                
                # Determine the amount of oil to pour.
                # It's the minimum of what's in the source and the remaining capacity of the destination.
                amount = min(next_state_list[i], capacities[j] - next_state_list[j])
                
                # If amount is 0, this pour does nothing, so skip it.
                if amount == 0:
                    continue
                
                # Perform the pour.
                next_state_list[i] -= amount
                next_state_list[j] += amount
                
                # Convert the new state back to a tuple to add it to the 'visited' set.
                next_state = tuple(next_state_list)
                
                if next_state not in visited:
                    visited.add(next_state)
                    # Create the string for this operation, e.g., "P(X, B)".
                    op_string = f"P({container_names[i]}, {container_names[j]})"
                    # Append this state and its path to the queue.
                    new_path = path + [op_string]
                    queue.append((next_state, new_path))
                    
    # If the queue becomes empty and the goal was not found.
    print("No solution found.")

if __name__ == "__main__":
    solve_oil_division()