from collections import deque

def solve_oil_puzzle():
    """
    Solves the oil pouring puzzle using Breadth-First Search (BFS)
    to find the shortest sequence of operations.
    """
    # Capacities of containers
    capacities = {'X': float('inf'), 'A': 8, 'B': 17, 'C': 21}
    cap_tuple = (float('inf'), 8, 17, 21)

    # Initial state: (oil_in_X, oil_in_A, oil_in_B, oil_in_C)
    initial_state = (39, 0, 0, 0)
    
    # Target condition: B=13, C=13, X+A=13. We search for the first state that meets this.
    # In this specific solution, the state is (5, 8, 13, 13)
    target_state = (5, 8, 13, 13)

    # Queue for BFS: stores tuples of (state, path_of_actions)
    queue = deque([(initial_state, [])])
    
    # Set to keep track of visited states to avoid cycles
    visited = {initial_state}
    
    # Map container names to state tuple indices
    state_map = {'X': 0, 'A': 1, 'B': 2, 'C': 3}
    
    # All possible pouring actions (Source, Destination)
    # Cannot pour back to X
    actions = [
        ('X', 'A'), ('X', 'B'), ('X', 'C'),
        ('A', 'B'), ('A', 'C'),
        ('B', 'A'), ('B', 'C'),
        ('C', 'A'), ('C', 'B')
    ]

    path = None
    while queue:
        current_state, current_path = queue.popleft()

        if current_state == target_state:
            path = current_path
            break

        # Try all possible pouring actions from the current state
        for src_name, dest_name in actions:
            src_idx = state_map[src_name]
            dest_idx = state_map[dest_name]

            src_amount = current_state[src_idx]
            dest_amount = current_state[dest_idx]
            dest_capacity = cap_tuple[dest_idx]

            # If source is empty or destination is full, this pour is not possible
            if src_amount == 0 or dest_amount == dest_capacity:
                continue

            # Calculate amount to pour (either empty source or fill destination)
            pour_amount = min(src_amount, dest_capacity - dest_amount)
            
            # Create the next state
            next_state_list = list(current_state)
            next_state_list[src_idx] -= pour_amount
            next_state_list[dest_idx] += pour_amount
            next_state = tuple(next_state_list)

            if next_state not in visited:
                visited.add(next_state)
                new_path = current_path + [f"P({src_name}, {dest_name})"]
                queue.append((next_state, new_path))
    
    # Print the solution
    if path:
        print("Solution found!")
        s = initial_state
        print(f"Start: X={s[0]}, A={s[1]}, B={s[2]}, C={s[3]}")
        for i, action in enumerate(path):
            src_name = action[2]
            dest_name = action[5]
            src_idx = state_map[src_name]
            dest_idx = state_map[dest_name]
            
            src_amount = s[src_idx]
            dest_amount = s[dest_idx]
            dest_capacity = cap_tuple[dest_idx]
            pour_amount = min(src_amount, dest_capacity - dest_amount)
            
            s_list = list(s)
            s_list[src_idx] -= pour_amount
            s_list[dest_idx] += pour_amount
            s = tuple(s_list)

            print(f"Step {i+1}: {action} -> X={s[0]}, A={s[1]}, B={s[2]}, C={s[3]}")
        print("\nFinal State Achieved:")
        print(f"Container B: {s[2]}L")
        print(f"Container C: {s[3]}L")
        print(f"Containers X+A: {s[0]}L + {s[1]}L = {s[0]+s[1]}L")
        print(f"\nThe last operation is: {path[-1]}")
    else:
        print("No solution found for the given target.")

solve_oil_puzzle()
<<<E>>>