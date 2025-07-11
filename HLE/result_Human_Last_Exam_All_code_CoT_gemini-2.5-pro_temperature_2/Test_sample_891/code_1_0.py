import collections

def solve_oil_division():
    """
    Solves the oil division puzzle using Breadth-First Search to find the shortest sequence.
    """
    # Define container capacities and total oil volume
    capacities = {'A': 8, 'B': 17, 'C': 21}
    container_names = ['A', 'B', 'C']
    initial_oil_in_X = 39

    # The initial state is (oil_in_A, oil_in_B, oil_in_C)
    initial_state = (0, 0, 0)

    # The queue will store tuples of (state, path_to_state)
    queue = collections.deque([(initial_state, [])])
    
    # The 'visited' set will store states we have already queued to avoid cycles
    visited = {initial_state}

    while queue:
        current_state, path = queue.popleft()
        a, b, c = current_state

        # --- GOAL CHECK ---
        # The goal is reached when B and C have 13L each.
        # The third 13L portion is then A + X.
        if b == 13 and c == 13:
            print("--- Solution Found ---")
            print(f"Shortest sequence has {len(path)} steps.")
            
            # Print the state transitions for clarity
            oil = {'X': initial_oil_in_X, 'A': 0, 'B': 0, 'C': 0}
            print(f"Start: X={oil['X']}, A={oil['A']}, B={oil['B']}, C={oil['C']}")
            
            for move in path:
                source, dest = move[2], move[4]
                if source == 'X':
                    poured = capacities[dest] - oil[dest]
                    oil['X'] -= poured
                    oil[dest] += poured
                else:
                    poured = min(oil[source], capacities[dest] - oil[dest])
                    oil[source] -= poured
                    oil[dest] += poured
                print(f"{move}: X={oil['X']}, A={oil['A']}, B={oil['B']}, C={oil['C']}")
                
            final_a, final_b, final_c = current_state
            final_x = initial_oil_in_X - (final_a + final_b + final_c)
            print("\n--- Final Distribution ---")
            print(f"Portion 1 (in Container B): {final_b}L")
            print(f"Portion 2 (in Container C): {final_c}L")
            print(f"Portion 3 (in A + X): {final_a}L + {final_x}L = {final_a + final_x}L")
            
            last_operation = path[-1]
            print(f"\nThe last operation is: {last_operation}")
            return

        # --- GENERATE NEXT STATES ---
        current_levels = {'A': a, 'B': b, 'C': c}
        oil_in_x = initial_oil_in_X - (a + b + c)
        
        # Action: Pour from X
        for dest in container_names:
            if current_levels[dest] < capacities[dest]:
                amount_to_fill = capacities[dest] - current_levels[dest]
                if oil_in_x >= amount_to_fill:
                    new_levels = current_levels.copy()
                    new_levels[dest] = capacities[dest]
                    new_state = (new_levels['A'], new_levels['B'], new_levels['C'])
                    if new_state not in visited:
                        visited.add(new_state)
                        queue.append((new_state, path + [f"P(X, {dest})"]))

        # Action: Pour between A, B, C
        for source in container_names:
            for dest in container_names:
                if source == dest:
                    continue
                if current_levels[source] > 0 and current_levels[dest] < capacities[dest]:
                    poured = min(current_levels[source], capacities[dest] - current_levels[dest])
                    new_levels = current_levels.copy()
                    new_levels[source] -= poured
                    new_levels[dest] += poured
                    new_state = (new_levels['A'], new_levels['B'], new_levels['C'])
                    if new_state not in visited:
                        visited.add(new_state)
                        queue.append((new_state, path + [f"P({source}, {dest})"]))

# Run the solver
solve_oil_division()