import collections

def solve_oil_division():
    """
    This function solves the oil division puzzle using Breadth-First Search (BFS)
    to find the shortest sequence of pouring actions.
    """
    # Capacities of containers A, B, C. X is considered to have infinite capacity.
    capacities = {'A': 8, 'B': 17, 'C': 21}
    # Initial state of oil in (X, A, B, C)
    start_state = (39, 0, 0, 0)
    # Target state where oil is divided into three 13L portions
    target_state = (13, 0, 13, 13)

    # The queue for BFS will store tuples of (state, path)
    # state: a tuple (oil_x, oil_a, oil_b, oil_c)
    # path: a list of action strings, e.g., ["P(X, B)", "P(B, C)"]
    queue = collections.deque([(start_state, [])])

    # A set to store visited states to prevent cycles and redundant explorations
    visited = {start_state}

    # Helper maps for container names and their corresponding index in the state tuple
    container_names = ['X', 'A', 'B', 'C']
    container_indices = {name: i for i, name in enumerate(container_names)}

    while queue:
        current_state, path = queue.popleft()

        # Check if the current state is the target state
        if current_state == target_state:
            print("Solution found!")
            print(f"The shortest sequence requires {len(path)} steps.")
            print("Sequence of operations:")
            print(" -> ".join(path))
            print("\nFinal distribution of oil:")
            
            x, a, b, c = current_state
            print(f"Container X: {x} liters")
            print(f"Container A: {a} liters")
            print(f"Container B: {b} liters")
            print(f"Container C: {c} liters")

            print("\nThe final equation is:")
            print(f"X ({x}) + A ({a}) + B ({b}) + C ({c}) = {x + a + b + c} liters")

            last_operation = path[-1]
            print(f"\nThe last operation to achieve the goal is: {last_operation}")

            # Map the operation string to the corresponding answer choice letter
            options = {
                'P(A, B)': 'A', 'P(A, C)': 'B', 'P(B, C)': 'C', 'P(C, B)': 'D',
                'P(C, A)': 'E', 'P(X, A)': 'F', 'P(X, B)': 'G', 'P(B, A)': 'H'
            }
            answer_letter = options.get(last_operation, "Unknown")
            print(f"<<<{answer_letter}>>>")
            return

        # Generate all valid next states by simulating all possible pours
        for src_name in container_names:
            # Destination cannot be the same as the source or container X
            for dest_name in ['A', 'B', 'C']:
                if src_name == dest_name:
                    continue

                src_idx = container_indices[src_name]
                dest_idx = container_indices[dest_name]

                current_src_amount = current_state[src_idx]
                current_dest_amount = current_state[dest_idx]

                # Cannot pour from an empty container
                if current_src_amount == 0:
                    continue
                
                dest_capacity = capacities[dest_name]
                # Cannot pour to a full container
                if current_dest_amount == dest_capacity:
                    continue

                # Calculate the amount to pour according to the rules
                pour_amount = min(current_src_amount, dest_capacity - current_dest_amount)

                # A pour must transfer a non-zero amount of oil
                if pour_amount == 0:
                    continue

                # Create the new state after the pour
                next_state_list = list(current_state)
                next_state_list[src_idx] -= pour_amount
                next_state_list[dest_idx] += pour_amount
                next_state = tuple(next_state_list)

                # If the new state has not been visited, add it to the queue
                if next_state not in visited:
                    visited.add(next_state)
                    new_path = path + [f"P({src_name}, {dest_name})"]
                    queue.append((next_state, new_path))
    
    print("No solution could be found.")

if __name__ == '__main__':
    solve_oil_division()