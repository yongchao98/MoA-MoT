from collections import deque

def solve_oil_puzzle():
    """
    Solves the oil pouring puzzle using Breadth-First Search (BFS) to find the
    shortest sequence of operations.
    """
    # Define container capacities
    capacities = {'A': 8, 'B': 17, 'C': 21}
    cap_tuple = (capacities['A'], capacities['B'], capacities['C'])
    
    # Initial state: (A=0, B=0, C=0), implies X=39
    initial_state_abc = (0, 0, 0)
    
    # Target state: (A=0, B=13, C=13), implies X=13
    target_state_abc = (0, 13, 13)

    # The queue will store tuples of: (current_state, path_to_current_state)
    queue = deque([(initial_state_abc, [])])
    
    # A set to store visited states to prevent cycles
    visited = {initial_state_abc}

    while queue:
        current_state, path = queue.popleft()
        a, b, c = current_state

        # Check if the current state is the target state
        if current_state == target_state_abc:
            # --- Solution Found ---
            print("The shortest sequence of operations is found:")
            
            # Print the step-by-step process
            temp_state = [39, 0, 0, 0] # X, A, B, C
            print(f"Initial State:\t(X={temp_state[0]}, A={temp_state[1]}, B={temp_state[2]}, C={temp_state[3]})")

            name_to_idx = {'X': 0, 'A': 1, 'B': 2, 'C': 3}
            full_caps = (float('inf'), capacities['A'], capacities['B'], capacities['C'])

            for i, action_str in enumerate(path):
                src_name = action_str[2]
                dst_name = action_str[4]
                src_idx = name_to_idx[src_name]
                dst_idx = name_to_idx[dst_name]

                # Calculate amount poured for this step
                if src_name == 'X':
                    # Pouring from X is a 'fill' operation
                    amount = full_caps[dst_idx] - temp_state[dst_idx]
                else:
                    # Pouring between A, B, C stops when source is empty or dest is full
                    amount = min(temp_state[src_idx], full_caps[dst_idx] - temp_state[dst_idx])
                
                temp_state[src_idx] -= amount
                temp_state[dst_idx] += amount
                
                print(f"Step {i+1}: {action_str}\t-> (X={temp_state[0]}, A={temp_state[1]}, B={temp_state[2]}, C={temp_state[3]})")
            
            print("\nFinal Answer:")
            print(f"The final operation to achieve the goal (X=13, A=0, B=13, C=13) is: {path[-1]}")
            return

        # --- Generate Next States ---
        x = 39 - (a + b + c)
        
        # Action Type 1: Pour from X to fill A, B, or C
        # P(X, A)
        if a < cap_tuple[0] and x >= (cap_tuple[0] - a):
            next_state = (cap_tuple[0], b, c)
            if next_state not in visited:
                visited.add(next_state)
                queue.append((next_state, path + ["P(X, A)"]))
        # P(X, B)
        if b < cap_tuple[1] and x >= (cap_tuple[1] - b):
            next_state = (a, cap_tuple[1], c)
            if next_state not in visited:
                visited.add(next_state)
                queue.append((next_state, path + ["P(X, B)"]))
        # P(X, C)
        if c < cap_tuple[2] and x >= (cap_tuple[2] - c):
            next_state = (a, b, cap_tuple[2])
            if next_state not in visited:
                visited.add(next_state)
                queue.append((next_state, path + ["P(X, C)"]))

        # Action Type 2: Pour between containers A, B, C
        containers = list(current_state)
        container_names = ['A', 'B', 'C']
        
        for i in range(3):  # Source container index
            for j in range(3):  # Destination container index
                if i == j:
                    continue
                
                # Check if a pour is possible
                if containers[i] > 0 and containers[j] < cap_tuple[j]:
                    amount = min(containers[i], cap_tuple[j] - containers[j])
                    
                    next_containers = list(containers)
                    next_containers[i] -= amount
                    next_containers[j] += amount
                    
                    next_state = tuple(next_containers)
                    if next_state not in visited:
                        visited.add(next_state)
                        action_str = f"P({container_names[i]}, {container_names[j]})"
                        queue.append((next_state, path + [action_str]))

    print("No solution found.")

# Execute the solver
solve_oil_puzzle()
<<<B>>>