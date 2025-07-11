import collections

def solve_oil_division():
    """
    Solves the oil division puzzle using Breadth-First Search (BFS)
    to find the shortest sequence of pours.
    """
    # Capacities of containers A, B, C. X is the source.
    capacities = {'X': float('inf'), 'A': 8, 'B': 17, 'C': 21}
    
    # Initial state: (oil_in_X, oil_in_A, oil_in_B, oil_in_C)
    initial_state = (39, 0, 0, 0)

    # The goal is to have three 13L portions.
    # This is achieved when container B and C have 13L each,
    # and the remaining oil (in X and A) also sums to 13L.
    
    # Queue for BFS: stores tuples of (state, path_list)
    queue = collections.deque([(initial_state, [])])
    
    # Set to keep track of visited states to avoid cycles
    visited = {initial_state}

    while queue:
        current_state, path = queue.popleft()
        
        x, a, b, c = current_state

        # Check if the current state meets the objective
        if b == 13 and c == 13:
            # The remaining oil is 39 - 13 - 13 = 13.
            # This state satisfies the condition of three 13L portions.
            print("Found the shortest solution!")
            print(f"The total number of steps is: {len(path)}")
            
            # Print the final state and the breakdown of the three portions
            print("\n--- Final State ---")
            print(f"Final distribution (X, A, B, C): {current_state}")
            print("The three 13L portions are:")
            print(f"1. Container B: {b}L")
            print(f"2. Container C: {c}L")
            # The final equation showing the third portion
            print(f"3. Container X + Container A: {x}L + {a}L = {x+a}L")
            
            last_operation = path[-1]
            print(f"\nThe last operation to achieve the goal is: {last_operation}")
            return

        # Generate next possible states by pouring
        containers = {'X': x, 'A': a, 'B': b, 'C': c}
        name_to_idx = {'X': 0, 'A': 1, 'B': 2, 'C': 3}
        
        # Source containers can be X, A, B, C
        for src_name in ['X', 'A', 'B', 'C']:
            src_amount = containers[src_name]
            if src_amount == 0:
                continue

            # Destination containers can be A, B, C
            for dest_name in ['A', 'B', 'C']:
                if src_name == dest_name:
                    continue

                dest_amount = containers[dest_name]
                dest_capacity = capacities[dest_name]
                space_in_dest = dest_capacity - dest_amount

                if space_in_dest == 0:
                    continue

                pour_amount = min(src_amount, space_in_dest)
                
                next_state_list = list(current_state)
                src_idx = name_to_idx[src_name]
                dest_idx = name_to_idx[dest_name]

                next_state_list[src_idx] -= pour_amount
                next_state_list[dest_idx] += pour_amount
                
                next_state = tuple(next_state_list)

                if next_state not in visited:
                    visited.add(next_state)
                    new_path = path + [f"P({src_name}, {dest_name})"]
                    queue.append((next_state, new_path))

    print("No solution found.")

# Execute the function to find the solution
solve_oil_division()