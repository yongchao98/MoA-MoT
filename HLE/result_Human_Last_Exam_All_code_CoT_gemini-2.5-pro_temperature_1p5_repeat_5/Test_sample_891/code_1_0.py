import collections

def solve_pouring_problem():
    """
    Solves the oil division puzzle using Breadth-First Search (BFS)
    to find the shortest sequence of pours.
    """
    
    # Define container capacities and names
    capacities = {'X': float('inf'), 'A': 8, 'B': 17, 'C': 21}
    container_names = ['X', 'A', 'B', 'C']
    
    # Define initial and goal states as tuples: (oil_X, oil_A, oil_B, oil_C)
    initial_state = (39, 0, 0, 0)
    goal_state = (13, 0, 13, 13)

    # Initialize the BFS queue with the starting state and an empty path
    # Each item in the queue is a tuple: (state, path_list)
    queue = collections.deque([(initial_state, [])])
    
    # Keep track of visited states to avoid cycles and redundant exploration
    visited = {initial_state}
    
    # Helper function to print the full solution path once found
    def print_solution(path):
        state = [39, 0, 0, 0]
        container_map_to_idx = {'X': 0, 'A': 1, 'B': 2, 'C': 3}
        
        print("Shortest sequence of pours found:")
        print(f"Start         -> (X: {state[0]:2d}, A: {state[1]:2d}, B: {state[2]:2d}, C: {state[3]:2d})")
        
        for i, action in enumerate(path):
            source_name = action[2]
            dest_name = action[5]
            
            s_idx = container_map_to_idx[source_name]
            d_idx = container_map_to_idx[dest_name]
            
            # Calculate amount to pour based on rules
            amount = min(state[s_idx], capacities[dest_name] - state[d_idx])
            
            # Update state
            state[s_idx] -= amount
            state[d_idx] += amount
            
            print(f"Step {i+1}: {action} -> (X: {state[0]:2d}, A: {state[1]:2d}, B: {state[2]:2d}, C: {state[3]:2d})")
    
    # Start the BFS loop
    while queue:
        current_state, path = queue.popleft()
        
        # Check if the goal state has been reached
        if current_state == goal_state:
            print_solution(path)
            last_operation = path[-1]
            
            print(f"\nThe goal is achieved. The final state is (X:13, A:0, B:13, C:13).")
            print(f"The last operation in the sequence is: {last_operation}")

            choices = {
                'A': 'P(A, B)', 'B': 'P(A, C)', 'C': 'P(B, C)', 'D': 'P(C, B)',
                'E': 'P(C, A)', 'F': 'P(X, A)', 'G': 'P(X, B)', 'H': 'P(B, A)'
            }
            
            final_answer_key = next((key for key, value in choices.items() if value == last_operation), None)
            
            print(f"This corresponds to answer choice: {final_answer_key}")
            print(f"<<<{final_answer_key}>>>")
            return

        # Generate all possible next states from the current state
        for s_idx, source_name in enumerate(container_names):
            for d_idx, dest_name in enumerate(container_names):
                # Rules: cannot pour into itself, cannot pour back to X
                if s_idx == d_idx or dest_name == 'X':
                    continue

                source_oil = current_state[s_idx]
                dest_oil = current_state[d_idx]
                dest_capacity = capacities[dest_name]

                # A pour is only possible if the source is not empty and the destination is not full
                if source_oil == 0 or dest_oil == dest_capacity:
                    continue
                
                # Calculate the amount to pour
                amount_to_pour = min(source_oil, dest_capacity - dest_oil)
                
                # Create the next state by applying the pour
                next_state_list = list(current_state)
                next_state_list[s_idx] -= amount_to_pour
                next_state_list[d_idx] += amount_to_pour
                next_state = tuple(next_state_list)
                
                # If the resulting state has not been visited, add it to the queue
                if next_state not in visited:
                    visited.add(next_state)
                    action_string = f"P({source_name}, {dest_name})"
                    new_path = path + [action_string]
                    queue.append((next_state, new_path))
    
    # If the queue becomes empty and the goal was not found
    print("No solution found.")

if __name__ == '__main__':
    solve_pouring_problem()