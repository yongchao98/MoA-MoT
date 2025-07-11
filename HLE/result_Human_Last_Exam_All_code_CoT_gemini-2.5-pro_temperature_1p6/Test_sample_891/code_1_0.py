from collections import deque

def solve_oil_division():
    """
    Solves the oil division puzzle using Breadth-First Search (BFS)
    to find the shortest sequence of pours.
    """
    
    # Capacities of containers A, B, C
    capacities = {'A': 8, 'B': 17, 'C': 21}
    initial_oil_x = 39

    # Target state for (A, B, C). X will be 13, and A+B+C = 26.
    # A cannot hold 13, so B and C must. A must be 0.
    target_state = (0, 13, 13)
    
    # Queue for BFS: stores ((state_a, state_b, state_c), path_list)
    # The state tracks A, B, C. X can be derived.
    queue = deque([((0, 0, 0), [])]) 
    
    # Set to keep track of visited states (a, b, c) to avoid cycles
    visited = set([(0, 0, 0)])

    while queue:
        (a, b, c), path = queue.popleft()

        # Check if we reached the target state
        if (a, b, c) == target_state:
            # We found the shortest path.
            print("Shortest sequence of operations found:")
            
            # To reconstruct the states for printing, we start from the beginning
            ca, cb, cc = 0, 0, 0
            cx = 39
            print(f"Start: X={cx:2d}, A={ca:2d}, B={cb:2d}, C={cc:2d}")
            
            for step in path:
                # This part is for detailed printing and not essential for the core BFS logic.
                s, d = step.replace("P(", "").replace(")", "").split(", ")
                
                # Update previous state to current state for the next iteration
                pa, pb, pc = ca, cb, cc
                px = cx
                
                if s == 'X':
                    src_val = px
                else:
                    src_val = {'A': pa, 'B': pb, 'C': pc}[s]

                dest_val = {'A': pa, 'B': pb, 'C': pc}[d]
                dest_cap = capacities[d]
                
                pour = min(src_val, dest_cap - dest_val)

                if s == 'X':
                    cx -= pour
                elif s == 'A':
                    ca -= pour
                elif s == 'B':
                    cb -= pour
                elif s == 'C':
                    cc -= pour

                if d == 'A':
                    ca += pour
                elif d == 'B':
                    cb += pour
                elif d == 'C':
                    cc += pour
                
                print(f"{len(path)-list(reversed(path)).index(step)}. {step:7s} -> X={cx:2d}, A={ca:2d}, B={cb:2d}, C={cc:2d}")

            print(f"\nFinal State: X={cx}L, B={cb}L, C={cc}L. (3 portions of 13L)")
            print(f"The last operation is: {path[-1]}")
            # The answer choices are A, B, C, etc. We map the last operation to the choice.
            options = {"P(A, B)": "A", "P(A, C)": "B", "P(B, C)": "C", "P(C, B)": "D", "P(C, A)": "E", "P(X, A)": "F", "P(X, B)": "G", "P(B, A)": "H"}
            answer = options.get(path[-1], "Unknown")
            print(f"\n<<<>>>") # Helper for final answer extraction.
            print(answer)
            return

        # --- Generate next possible states ---
        current_x = initial_oil_x - (a + b + c)
        
        # All possible pours: from a source 's' to a destination 'd'
        # Sources can be X, A, B, C. Destinations can be A, B, C.
        sources = ['X', 'A', 'B', 'C']
        destinations = ['A', 'B', 'C']

        for s_name in sources:
            for d_name in destinations:
                if s_name == d_name:
                    continue

                # Get current values
                state_tuple = (a, b, c)
                
                if s_name == 'X':
                    s_val = current_x
                else:
                    s_val = state_tuple[ord(s_name) - ord('A')]

                if s_val == 0:
                    continue

                d_idx = ord(d_name) - ord('A')
                d_val = state_tuple[d_idx]
                d_cap = capacities[d_name]
                
                if d_val == d_cap:
                    continue
                
                # Perform the pour
                pour_amount = min(s_val, d_cap - d_val)
                
                new_state_list = list(state_tuple)
                
                if s_name != 'X':
                    s_idx = ord(s_name) - ord('A')
                    new_state_list[s_idx] -= pour_amount
                
                new_state_list[d_idx] += pour_amount
                new_state = tuple(new_state_list)

                # Constraint: X must not drop below 13
                new_x = initial_oil_x - sum(new_state)
                if new_x < 13:
                    continue

                if new_state not in visited:
                    visited.add(new_state)
                    new_path = path + [f"P({s_name}, {d_name})"]
                    queue.append((new_state, new_path))
                    
    print("No solution found.")

solve_oil_division()