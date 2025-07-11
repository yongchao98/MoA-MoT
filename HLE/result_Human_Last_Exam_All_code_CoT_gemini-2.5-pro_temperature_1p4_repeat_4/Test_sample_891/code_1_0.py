import collections

def solve_oil_division_puzzle():
    """
    Finds the shortest sequence of pouring actions to divide 39L of oil
    into three 13L portions using Breadth-First Search (BFS).
    """
    # Capacities of containers (X, A, B, C)
    capacities = (float('inf'), 8, 17, 21)
    # Container names for printing the path
    container_names = ['X', 'A', 'B', 'C']

    # Initial state: (oil_in_X, oil_in_A, oil_in_B, oil_in_C)
    initial_state = (39, 0, 0, 0)
    # Goal state: three portions of 13L, which must be (13, 0, 13, 13)
    target_state = (13, 0, 13, 13)

    # The queue for BFS will store tuples of (state, path_to_state)
    queue = collections.deque([(initial_state, [])])
    # The 'visited' set stores states we've already processed to avoid cycles
    visited = {initial_state}

    while queue:
        current_state, path = queue.popleft()

        # Check if we have reached the goal state
        if current_state == target_state:
            # Found the shortest path. Get the last operation.
            last_operation = path[-1]
            
            # To fulfill the "output each number in the final equation" requirement,
            # we determine the state just before the final move (penultimate state).
            penultimate_state = initial_state
            for move in path[:-1]:
                source_name = move[2]
                dest_name = move[5]
                s_idx = container_names.index(source_name)
                d_idx = container_names.index(dest_name)
                
                s_vol = penultimate_state[s_idx]
                d_vol = penultimate_state[d_idx]
                d_cap = capacities[d_idx]
                
                pour = min(s_vol, d_cap - d_vol)
                
                next_vols = list(penultimate_state)
                next_vols[s_idx] -= pour
                next_vols[d_idx] += pour
                penultimate_state = tuple(next_vols)

            print(f"The shortest sequence has {len(path)} steps.")
            print(f"The last operation is: {last_operation}")
            print(f"This operation transforms the state from:")
            print(f"X={penultimate_state[0]}, A={penultimate_state[1]}, B={penultimate_state[2]}, C={penultimate_state[3]}")
            print("to the final state:")
            print(f"X={target_state[0]}, A={target_state[1]}, B={target_state[2]}, C={target_state[3]}")
            print("\nThe numbers in this final state transition are:")
            print(f"{penultimate_state[0]}, {penultimate_state[1]}, {penultimate_state[2]}, {penultimate_state[3]} -> {target_state[0]}, {target_state[1]}, {target_state[2]}, {target_state[3]}")
            return

        # If not the goal, generate all possible next states
        for s_idx in range(4):  # Source container index
            # Destination cannot be X (no pouring back to X)
            for d_idx in range(1, 4):
                if s_idx == d_idx:
                    continue

                s_vol = current_state[s_idx]
                d_vol = current_state[d_idx]
                d_cap = capacities[d_idx]

                # Pouring is only possible if source is not empty and destination is not full
                if s_vol > 0 and d_vol < d_cap:
                    pour_amount = min(s_vol, d_cap - d_vol)
                    
                    next_vols = list(current_state)
                    next_vols[s_idx] -= pour_amount
                    next_vols[d_idx] += pour_amount
                    next_state = tuple(next_vols)

                    if next_state not in visited:
                        visited.add(next_state)
                        move_str = f"P({container_names[s_idx]}, {container_names[d_idx]})"
                        new_path = path + [move_str]
                        queue.append((next_state, new_path))
    
    print("No solution was found.")

solve_oil_division_puzzle()