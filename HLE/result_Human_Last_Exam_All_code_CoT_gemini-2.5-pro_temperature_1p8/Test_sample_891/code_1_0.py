import collections

def solve_oil_puzzle():
    """
    Solves the oil division puzzle using Breadth-First Search (BFS)
    to find the shortest sequence of pours.
    """
    
    # Capacities of containers (A, B, C)
    capacities = (8, 17, 21)
    # Initial total oil
    initial_oil = 39

    # The goal state is (a=0, b=13, c=13)
    goal_state = (0, 13, 13)

    # Initial state (a=0, b=0, c=0) with an empty path
    initial_state = (0, 0, 0)
    
    # The queue for BFS will store tuples of (state, path)
    # where path is a list of action descriptions.
    queue = collections.deque([(initial_state, [])])
    
    # A set to keep track of visited states to avoid cycles and redundant computations.
    visited = {initial_state}

    while queue:
        current_state, path = queue.popleft()
        a, b, c = current_state

        # Check if the goal state is reached
        if current_state == goal_state:
            print("Shortest sequence found:")
            # Initial state printout
            initial_a, initial_b, initial_c = (0, 0, 0)
            initial_x = initial_oil - (initial_a + initial_b + initial_c)
            print(f"Start State: (X={initial_x}, A={initial_a}, B={initial_b}, C={initial_c})")

            # Print the sequence of operations
            for move in path:
                action, src, dest, state_after = move
                x_after = initial_oil - sum(state_after)
                print(f"{action}({src}, {dest}) -> State: (X={x_after}, A={state_after[0]}, B={state_after[1]}, C={state_after[2]})")
            
            # Extract and print the last operation for the final answer
            last_op, last_src, last_dest, _ = path[-1]
            print(f"\nThe final equation to reach the goal is by pouring from container {last_src} to {last_dest}.")
            print(f"The last operation is: {last_op}({last_src}, {last_dest})")
            return

        # ---- Generate all possible next states ----
        
        # Current amount in source container X
        x = initial_oil - (a + b + c)
        
        # State possibilities: tuple of (current_vols, source_name)
        states_to_pour_from = [
            ((a, b, c), "A", a),
            ((a, b, c), "B", b),
            ((a, b, c), "C", c),
        ]
        if x > 0:
            states_to_pour_from.append(((a,b,c), "X", x))

        # Destination possibilities: tuple of (dest_idx, dest_name, dest_vol)
        destinations = [(0, "A", a), (1, "B", b), (2, "C", c)]

        for state_tuple, src_name, src_vol in states_to_pour_from:
            if src_vol == 0:
                continue
            
            for dest_idx, dest_name, dest_vol in destinations:
                # Cannot pour into the same container or back to X
                if src_name == dest_name or dest_name == 'X':
                    continue

                dest_cap = capacities[dest_idx]
                
                # If destination is full, skip
                if dest_vol == dest_cap:
                    continue
                
                # Calculate amount to pour
                amount = min(src_vol, dest_cap - dest_vol)
                
                # Create the new state by transferring the oil
                next_state_list = list(current_state)
                
                # Update destination
                next_state_list[dest_idx] += amount
                
                # Update source
                if src_name == "A":
                    next_state_list[0] -= amount
                elif src_name == "B":
                    next_state_list[1] -= amount
                elif src_name == "C":
                    next_state_list[2] -= amount
                # Source X is handled by recalculating `x` in the next loop

                next_state = tuple(next_state_list)
                if next_state not in visited:
                    visited.add(next_state)
                    new_path = path + [(f"P", src_name, dest_name, next_state)]
                    queue.append((next_state, new_path))
    
    print("No solution found.")

solve_oil_puzzle()
<<<A>>>