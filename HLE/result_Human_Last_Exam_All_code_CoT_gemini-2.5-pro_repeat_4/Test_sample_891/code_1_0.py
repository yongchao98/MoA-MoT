from collections import deque

def solve_pouring_problem():
    """
    Solves the oil pouring problem using Breadth-First Search (BFS)
    to find the shortest sequence of operations.
    """
    # Capacities of the containers
    capacities = {'X': float('inf'), 'A': 8, 'B': 17, 'C': 21}
    
    # Initial state: (oil_in_X, oil_in_A, oil_in_B, oil_in_C)
    initial_state = (39, 0, 0, 0)
    
    # Target state
    target_state = (13, 0, 13, 13)

    # Queue for BFS: stores tuples of (state, path_list)
    queue = deque([(initial_state, [])])
    
    # A set to keep track of visited states to avoid cycles
    visited = {initial_state}

    while queue:
        current_state, path = queue.popleft()
        x, a, b, c = current_state

        # Check if we have reached the target state
        if current_state == target_state:
            # Print the final operation
            print(f"The last operation is: {path[-1]}")
            # Print the final state equation as requested
            print("Final container distribution:")
            print(f"Container X: {x}L")
            print(f"Container A: {a}L")
            print(f"Container B: {b}L")
            print(f"Container C: {c}L")
            print("Final equation:")
            print(f"{x} + {a} + {b} + {c} = {x + a + b + c}")
            return path

        # --- Generate all possible next states ---
        
        # Pour from X to A, B, C (Fill destination)
        # Note: The problem implies filling from X is a primary action.
        # We assume it's allowed even if the destination is not empty (topping up).
        sources = {'X': x, 'A': a, 'B': b, 'C': c}
        
        # Iterate through all possible (source, destination) pairs
        pour_pairs = [('X', 'A'), ('X', 'B'), ('X', 'C'),
                      ('A', 'B'), ('A', 'C'),
                      ('B', 'A'), ('B', 'C'),
                      ('C', 'A'), ('C', 'B')]

        for source_name, dest_name in pour_pairs:
            source_val = sources[source_name]
            dest_val = sources[dest_name]

            if source_val == 0:
                continue

            if source_name == 'X':
                # Pouring from X means filling the destination
                amount = capacities[dest_name] - dest_val
                if amount == 0: continue # Destination is already full
                if x < amount: continue # Not enough oil in X
                
                new_state_list = [x, a, b, c]
                new_state_list[0] -= amount
                idx_map = {'A':1, 'B':2, 'C':3}
                new_state_list[idx_map[dest_name]] += amount
                next_state = tuple(new_state_list)
            else:
                # Pouring between A, B, C
                dest_cap = capacities[dest_name]
                if dest_val == dest_cap: continue # Destination is already full

                amount = min(source_val, dest_cap - dest_val)
                
                new_state_list = [x, a, b, c]
                idx_map = {'A':1, 'B':2, 'C':3}
                new_state_list[idx_map[source_name]] -= amount
                new_state_list[idx_map[dest_name]] += amount
                next_state = tuple(new_state_list)

            if next_state not in visited:
                visited.add(next_state)
                new_path = path + [f"P({source_name}, {dest_name})"]
                queue.append((next_state, new_path))
                
    print("No solution found.")
    return None

# Execute the solver
solve_pouring_problem()