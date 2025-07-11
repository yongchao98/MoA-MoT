from collections import deque

def solve_oil_puzzle():
    """
    Solves the oil pouring puzzle using Breadth-First Search (BFS)
    to find the shortest sequence of operations.
    """
    # Capacities of containers (X is infinite, A, B, C)
    capacities = {'X': float('inf'), 'A': 8, 'B': 17, 'C': 21}
    container_names = ['X', 'A', 'B', 'C']
    
    # Map names to indices for state tuple
    name_to_idx = {name: i for i, name in enumerate(container_names)}
    idx_to_name = {i: name for i, name in enumerate(container_names)}

    # Initial state: (oil_in_X, oil_in_A, oil_in_B, oil_in_C)
    start_state = (39, 0, 0, 0)
    # Goal state
    target_state = (13, 0, 13, 13)

    # BFS queue stores (state, path_of_actions)
    queue = deque([(start_state, [])])
    # Visited set to store states we've already seen
    visited = {start_state}

    while queue:
        current_state, path = queue.popleft()

        # If we've reached the target, we're done
        if current_state == target_state:
            print("Shortest sequence found:")
            # We need to print each number in the final equation.
            # The final equation is the last pour.
            # Let's find the state before the last pour.
            
            # To reconstruct the state before the last pour, we can "un-pour".
            # However, it's easier to just trace the path forward to get the last state's values.
            
            # Let's trace the path to get the numbers for the last step
            state_trace = start_state
            for i, action in enumerate(path):
                source_name, dest_name = action[2], action[5]
                s_idx, d_idx = name_to_idx[source_name], name_to_idx[dest_name]
                
                s_cap_val = capacities[source_name] if source_name != 'X' else float('inf')
                d_cap_val = capacities[dest_name]

                amount_to_pour = min(state_trace[s_idx], d_cap_val - state_trace[d_idx])
                
                # If this is the last action, print the details
                if i == len(path) - 1:
                    print(f"Final Step: Pour {int(amount_to_pour)} liters from Container {source_name} (contains {int(state_trace[s_idx])}) to Container {dest_name} (contains {int(state_trace[d_idx])}).")
                    print(f"Equation: Container {dest_name} new amount = {int(state_trace[d_idx])} + {int(amount_to_pour)} = {int(state_trace[d_idx] + amount_to_pour)}")
                    print(f"Container {source_name} new amount = {int(state_trace[s_idx])} - {int(amount_to_pour)} = {int(state_trace[s_idx] - amount_to_pour)}")

                new_state_list = list(state_trace)
                new_state_list[s_idx] -= amount_to_pour
                new_state_list[d_idx] += amount_to_pour
                state_trace = tuple(new_state_list)

            print("\nFull Path:")
            for i, step in enumerate(path, 1):
                print(f"Step {i}: {step}")

            print(f"\nFinal State: (X:{state_trace[0]}, A:{state_trace[1]}, B:{state_trace[2]}, C:{state_trace[3]})")
            print(f"\nThe last operation is: {path[-1]}")
            return path[-1]

        # Generate next possible states from all valid pours
        for i in range(4):  # source index
            for j in range(4):  # destination index
                if i == j:
                    continue

                # Rule: Cannot pour back into X
                if j == 0:
                    continue

                source_name = idx_to_name[i]
                dest_name = idx_to_name[j]
                
                source_val = current_state[i]
                dest_val = current_state[j]
                dest_cap = capacities[dest_name]

                if source_val == 0 or dest_val == dest_cap:
                    continue

                amount = min(source_val, dest_cap - dest_val)

                if amount == 0:
                    continue
                
                next_state_list = list(current_state)
                next_state_list[i] -= amount
                next_state_list[j] += amount
                next_state = tuple(next_state_list)

                if next_state not in visited:
                    visited.add(next_state)
                    new_path = path + [f"P({source_name}, {dest_name})"]
                    queue.append((next_state, new_path))
    
    print("No solution found.")
    return None

solve_oil_puzzle()
<<<A>>>