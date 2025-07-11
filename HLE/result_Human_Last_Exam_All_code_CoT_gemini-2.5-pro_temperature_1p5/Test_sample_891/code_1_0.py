import collections

def solve_oil_division_problem():
    """
    This function solves the oil division puzzle by finding the shortest
    sequence of pouring actions using a Breadth-First Search (BFS) algorithm.
    """
    # Define container capacities
    capacities = {'A': 8, 'B': 17, 'C': 21}
    
    # The initial state of the containers (X, A, B, C)
    initial_state = (39, 0, 0, 0)

    # --- BFS setup ---
    # The queue will store tuples of (state, path)
    queue = collections.deque([(initial_state, [])])
    # The visited set stores states to prevent cycles
    visited = {initial_state}

    print("Searching for the shortest sequence of pouring actions...")

    # Start BFS loop
    while queue:
        current_state, path = queue.popleft()
        x, a, b, c = current_state

        # Check if the current state meets the target condition
        if b == 13 and c == 13 and (x + a) == 13:
            print("\n--- Solution Found! ---")
            print(f"The shortest path has {len(path)} steps.")
            
            # Print the entire sequence for clarity
            temp_state = list(initial_state)
            name_to_idx = {'X': 0, 'A': 1, 'B': 2, 'C': 3}
            print(f"Step 0: Start          -> (X:{temp_state[0]}, A:{temp_state[1]}, B:{temp_state[2]}, C:{temp_state[3]})")
            
            for i, move in enumerate(path):
                s_name = move[2]
                d_name = move[5]
                s_idx = name_to_idx[s_name]
                d_idx = name_to_idx[d_name]
                d_cap = capacities[d_name]

                pour = min(temp_state[s_idx], d_cap - temp_state[d_idx])
                temp_state[s_idx] -= pour
                temp_state[d_idx] += pour
                print(f"Step {i+1}: {move.ljust(10)} -> (X:{temp_state[0]}, A:{temp_state[1]}, B:{temp_state[2]}, C:{temp_state[3]})")
            
            print("\n--- Final Answer ---")
            last_operation = path[-1]
            # This part prints the information related to the last step, as requested.
            print(f"The state before the last operation was (X:{x}, A:{a}, B:{b}, C:{c}).")
            print(f"The last operation is: {last_operation}")
            return

        # --- Generate next possible states ---
        # Map container names to their current volumes for easier access
        vols = {'X': x, 'A': a, 'B': b, 'C': c}
        sources = ['X', 'A', 'B', 'C']
        destinations = ['A', 'B', 'C']

        for s_name in sources:
            for d_name in destinations:
                if s_name == d_name:
                    continue

                s_vol = vols[s_name]
                d_vol = vols[d_name]
                d_cap = capacities[d_name]

                # A pour is only possible if the source is not empty and destination is not full
                if s_vol > 0 and d_vol < d_cap:
                    pour_amount = min(s_vol, d_cap - d_vol)
                    
                    next_vols = vols.copy()
                    next_vols[s_name] -= pour_amount
                    next_vols[d_name] += pour_amount
                    
                    next_state = (next_vols['X'], next_vols['A'], next_vols['B'], next_vols['C'])

                    if next_state not in visited:
                        visited.add(next_state)
                        new_path = path + [f"P({s_name}, {d_name})"]
                        queue.append((next_state, new_path))
    
    print("No solution was found.")

solve_oil_division_problem()
<<<H>>>