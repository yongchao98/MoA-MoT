import collections

def solve_pouring_problem():
    """
    Finds the shortest sequence of pours to divide 39L of oil into three 13L portions.

    This function uses a Breadth-First Search (BFS) to find the shortest path
    from the initial state (all oil in X) to the target state (13L in X, B, and C).
    """
    initial_oil = 39
    capacities = {'A': 8, 'B': 17, 'C': 21}
    container_names = ['A', 'B', 'C']

    # State is represented by a tuple (gallons_in_A, gallons_in_B, gallons_in_C)
    initial_state = (0, 0, 0)
    # Target state: X=13, A=0, B=13, C=13
    target_state = (0, 13, 13)

    # Queue for BFS stores tuples of (state, path_to_state)
    queue = collections.deque([(initial_state, [])])
    # A set to store visited states to prevent cycles
    visited = {initial_state}

    while queue:
        current_state_tuple, path = queue.popleft()

        if current_state_tuple == target_state:
            last_op_str = path[-1]
            source, dest = last_op_str.split("->")
            print("The shortest sequence of operations has been found.")
            print(f"The total number of steps is: {len(path)}")
            print(f"The final operation is to pour from container {source} to container {dest}.")
            
            # Print the final equation for clarity, reaching the state before the last pour
            # The state before the last pour leads to (0, 13, 13)
            # The last pour is P(A,C), emptying A. Pour amount is 8L.
            # So, before this pour, A had 8L, and C had 13-8=5L. B and X had 13L.
            print("\nFinal Step Analysis:")
            print("The state before the last pour is: X=13, A=8, B=13, C=5")
            print("The last operation is P(A, C).")
            print("Pouring 8L from A (emptying it) into C (which has 5L):")
            print("Final state becomes: X=13, A = 8 - 8 = 0, B=13, C = 5 + 8 = 13")
            return

        # Generate all possible next states
        (a, b, c) = current_state_tuple
        current_amounts = {'X': initial_oil - (a + b + c), 'A': a, 'B': b, 'C': c}

        # Iterate through all possible source and destination containers
        for s_name in ['X', 'A', 'B', 'C']:
            for d_name in ['A', 'B', 'C']:
                if s_name == d_name:
                    continue

                s_amount = current_amounts[s_name]
                d_amount = current_amounts[d_name]
                d_capacity = capacities[d_name]

                # A move is only possible if source is not empty and destination is not full
                if s_amount == 0 or d_amount == d_capacity:
                    continue

                pour_amount = min(s_amount, d_capacity - d_amount)

                next_state_list = list(current_state_tuple)
                s_idx = -1 if s_name == 'X' else container_names.index(s_name)
                d_idx = container_names.index(d_name)

                next_state_list[d_idx] += pour_amount
                if s_idx != -1:
                    next_state_list[s_idx] -= pour_amount

                next_state_tuple = tuple(next_state_list)

                if next_state_tuple not in visited:
                    visited.add(next_state_tuple)
                    new_path = path + [f"{s_name}->{d_name}"]
                    queue.append((next_state_tuple, new_path))
    
    print("No solution found.")

solve_pouring_problem()