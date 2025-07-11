import collections

def solve_oil_puzzle():
    """
    Solves the oil division puzzle using Breadth-First Search (BFS)
    to find the shortest sequence of pouring actions.
    """
    # Capacities of the containers
    cap_A = 8
    cap_B = 17
    cap_C = 21

    # Initial state (X, A, B, C) and target state
    initial_state = (39, 0, 0, 0)
    target_state = (13, 0, 13, 13)

    # Queue for BFS: stores tuples of (state, path_of_actions)
    queue = collections.deque([(initial_state, [])])
    
    # A set to store visited states to prevent cycles and redundant work
    visited = {initial_state}

    # Action to Choice mapping
    action_map = {
        "P(A, B)": "A", "P(A, C)": "B", "P(B, C)": "C",
        "P(C, B)": "D", "P(C, A)": "E", "P(X, A)": "F",
        "P(X, B)": "G", "P(B, A)": "H"
    }

    while queue:
        current_state, path = queue.popleft()
        x, a, b, c = current_state

        # Check if the target state has been reached
        if current_state == target_state:
            print("Shortest solution found in {} steps:".format(len(path)))
            
            # Print the sequence of states (the "final equation")
            temp_state = list(initial_state)
            print("Initial State : (X: {}, A: {}, B: {}, C: {})".format(temp_state[0], temp_state[1], temp_state[2], temp_state[3]))
            
            for i, action in enumerate(path):
                source_name, dest_name = action[2], action[5]
                
                containers = {'X': temp_state[0], 'A': temp_state[1], 'B': temp_state[2], 'C': temp_state[3]}
                caps = {'A': cap_A, 'B': cap_B, 'C': cap_C}

                source_val = containers[source_name]
                dest_val = containers[dest_name]
                dest_cap = caps[dest_name]

                pour_amount = min(source_val, dest_cap - dest_val)

                containers[source_name] -= pour_amount
                containers[dest_name] += pour_amount
                
                temp_state = [containers['X'], containers['A'], containers['B'], containers['C']]
                
                print("Step {:<2}: {} -> State: (X: {}, A: {}, B: {}, C: {})".format(
                    i + 1, action, temp_state[0], temp_state[1], temp_state[2], temp_state[3]
                ))

            print("\nThe final state is X=13, B=13, C=13 as required.")
            last_operation = path[-1]
            print(f"\nThe last operation is: {last_operation}")
            
            final_answer_choice = action_map.get(last_operation, "Unknown")
            print(f"<<<{final_answer_choice}>>>")
            return

        # --- Generate all possible next states ---
        # Note: Pouring from X is only allowed to fill the destination.
        # This is handled by the min() logic where source (X) is large enough.
        
        # Define all possible moves
        moves = []
        # Pour from X
        if x > 0:
            if a < cap_A: moves.append(('X', 'A', (x, a, b, c)))
            if b < cap_B: moves.append(('X', 'B', (x, a, b, c)))
            if c < cap_C: moves.append(('X', 'C', (x, a, b, c)))
        # Pour from A
        if a > 0:
            if b < cap_B: moves.append(('A', 'B', (x, a, b, c)))
            if c < cap_C: moves.append(('A', 'C', (x, a, b, c)))
        # Pour from B
        if b > 0:
            if a < cap_A: moves.append(('B', 'A', (x, a, b, c)))
            if c < cap_C: moves.append(('B', 'C', (x, a, b, c)))
        # Pour from C
        if c > 0:
            if a < cap_A: moves.append(('C', 'A', (x, a, b, c)))
            if b < cap_B: moves.append(('C', 'B', (x, a, b, c)))
        
        container_map = {'X': 0, 'A': 1, 'B': 2, 'C': 3}
        caps = {'A': cap_A, 'B': cap_B, 'C': cap_C}
        
        for source_name, dest_name, state in moves:
            s_idx, d_idx = container_map[source_name], container_map[dest_name]
            
            source_val = state[s_idx]
            dest_val = state[d_idx]
            dest_cap = caps[dest_name]

            # Calculate amount to pour
            pour = min(source_val, dest_cap - dest_val)

            if pour > 0:
                next_state_list = list(state)
                next_state_list[s_idx] -= pour
                next_state_list[d_idx] += pour
                next_state = tuple(next_state_list)

                if next_state not in visited:
                    visited.add(next_state)
                    new_path = path + ["P({}, {})".format(source_name, dest_name)]
                    queue.append((next_state, new_path))

    print("No solution found.")

if __name__ == "__main__":
    solve_oil_puzzle()