from collections import deque

def solve_oil_puzzle():
    """
    Solves the oil pouring puzzle using Breadth-First Search to find the shortest sequence of operations.
    """
    # Capacities of containers A, B, C
    capacities = {'A': 8, 'B': 17, 'C': 21}
    containers = ['A', 'B', 'C']
    
    # Initial state: (oil_a, oil_b, oil_c)
    initial_state = (0, 0, 0)
    # Goal state: (oil_a, oil_b, oil_c) corresponding to X=13, A=0, B=13, C=13
    goal_state = (0, 13, 13)
    total_oil = 39

    # Queue for BFS: stores tuples of (state, path_list)
    queue = deque([(initial_state, [])])
    # Set to keep track of visited states to avoid cycles
    visited = {initial_state}

    while queue:
        current_state, path = queue.popleft()

        # If the goal is found, print the last operation and terminate.
        if current_state == goal_state:
            print("Shortest sequence found:")
            # To match the answer choices format, we print just the last operation.
            last_operation = path[-1]
            print(f"The last operation is: {last_operation}")
            return

        # Unpack current amounts
        a, b, c = current_state
        current_levels = {'A': a, 'B': b, 'C': c}
        x = total_oil - (a + b + c)

        # --- Define all possible next states ---

        # Action Type 1: Pour from X to A, B, or C
        for dest_name in containers:
            dest_level = current_levels[dest_name]
            dest_cap = capacities[dest_name]
            
            if x > 0 and dest_level < dest_cap:
                pour_amount = min(x, dest_cap - dest_level)
                
                new_levels = current_levels.copy()
                new_levels[dest_name] += pour_amount
                new_state = (new_levels['A'], new_levels['B'], new_levels['C'])

                if new_state not in visited:
                    visited.add(new_state)
                    new_path = path + [f"P(X, {dest_name})"]
                    queue.append((new_state, new_path))
        
        # Action Type 2: Pour between A, B, and C
        for src_name in containers:
            for dest_name in containers:
                if src_name == dest_name:
                    continue

                src_level = current_levels[src_name]
                dest_level = current_levels[dest_name]
                dest_cap = capacities[dest_name]

                # Cannot pour from an empty container or into a full one
                if src_level == 0 or dest_level == dest_cap:
                    continue
                
                pour_amount = min(src_level, dest_cap - dest_level)
                
                new_levels = current_levels.copy()
                new_levels[src_name] -= pour_amount
                new_levels[dest_name] += pour_amount
                new_state = (new_levels['A'], new_levels['B'], new_levels['C'])

                if new_state not in visited:
                    visited.add(new_state)
                    new_path = path + [f"P({src_name}, {dest_name})"]
                    queue.append((new_state, new_path))

    print("No solution found.")

# Execute the solver
solve_oil_puzzle()