from collections import deque

def solve_oil_puzzle():
    """
    Solves the oil division puzzle using Breadth-First Search (BFS)
    to find the shortest sequence of pouring actions.
    """
    # Capacities of containers A, B, C
    caps = {'A': 8, 'B': 17, 'C': 21}
    container_names = ['A', 'B', 'C']
    container_caps = [caps[name] for name in container_names]

    # Initial state (oil in A, B, C)
    initial_state = (0, 0, 0)
    # Target state: A=0, B=13, C=13. (X will have 39 - 26 = 13)
    target_state = (0, 13, 13)

    # Queue for BFS: stores (state, path_list)
    queue = deque([(initial_state, [])])
    # Set to keep track of visited states to avoid cycles
    visited = {initial_state}

    while queue:
        current_state, path = queue.popleft()
        
        # Check if we reached the target state
        if current_state == target_state:
            # Reconstruct the state at each step for a clear trace
            print("Found a shortest solution. The sequence of operations is:")
            
            state_trace = []
            s = list(initial_state)
            state_trace.append( (39 - sum(s), s[0], s[1], s[2]) )

            for action in path:
                source, dest = action.replace("P(","").replace(")","").split(", ")
                
                temp_state = list(state_trace[-1][1:]) # Get last (a,b,c)
                
                if source == 'X':
                    dest_idx = container_names.index(dest)
                    pour_amount = container_caps[dest_idx] - temp_state[dest_idx]
                    temp_state[dest_idx] += pour_amount
                else:
                    src_idx = container_names.index(source)
                    dest_idx = container_names.index(dest)
                    
                    vol_src = temp_state[src_idx]
                    vol_dest = temp_state[dest_idx]
                    cap_dest = container_caps[dest_idx]
                    
                    if vol_src + vol_dest <= cap_dest: # Empty source
                        pour_amount = vol_src
                        temp_state[dest_idx] += pour_amount
                        temp_state[src_idx] = 0
                    else: # Fill destination
                        pour_amount = cap_dest - vol_dest
                        temp_state[dest_idx] += pour_amount
                        temp_state[src_idx] -= pour_amount
                
                s = tuple(temp_state)
                state_trace.append( (39 - sum(s), s[0], s[1], s[2]) )

            print(f"Initial State  : X=39, A=0, B=0, C=0")
            for i, p in enumerate(path):
                s_next = state_trace[i+1]
                print(f"{i+1}. {p:7s} -> X={s_next[0]:2d}, A={s_next[1]:2d}, B={s_next[2]:2d}, C={s_next[3]:2d}")

            last_op = path[-1]
            print(f"\nThe final state (X=13, A=0, B=13, C=13) is achieved.")
            print(f"The last operation is: {last_op}")
            
            choices = {
                "P(A, B)": "A", "P(A, C)": "B", "P(B, C)": "C", "P(C, B)": "D",
                "P(C, A)": "E", "P(X, A)": "F", "P(X, B)": "G", "P(B, A)": "H"
            }
            if last_op in choices:
                print(f"<<<{choices[last_op]}>>>")
            else:
                 print("Error: Last operation not found in choices.")
            return

        a, b, c = current_state
        x = 39 - (a + b + c)
        
        # --- Generate all possible next states ---

        # 1. Pour from X (only Fill Destination is possible)
        for i in range(3):
            name_dest = container_names[i]
            vol_dest = current_state[i]
            cap_dest = container_caps[i]
            
            pour_amount = cap_dest - vol_dest
            if x >= pour_amount and pour_amount > 0:
                next_state_list = list(current_state)
                next_state_list[i] += pour_amount
                next_state = tuple(next_state_list)
                if next_state not in visited:
                    visited.add(next_state)
                    queue.append((next_state, path + [f"P(X, {name_dest})"]))
        
        # 2. Pour between A, B, C
        for i in range(3):  # source index
            for j in range(3):  # destination index
                if i == j: continue

                vol_src = current_state[i]
                if vol_src == 0: continue
                
                name_src = container_names[i]
                name_dest = container_names[j]
                vol_dest = current_state[j]
                cap_dest = container_caps[j]

                # Action: Empty source into destination
                if vol_src + vol_dest <= cap_dest:
                    next_state_list = list(current_state)
                    next_state_list[j] += vol_src
                    next_state_list[i] = 0
                    next_state = tuple(next_state_list)
                    if next_state not in visited:
                        visited.add(next_state)
                        queue.append((next_state, path + [f"P({name_src}, {name_dest})"]))
                # Action: Fill destination from source
                else:
                    pour_amount = cap_dest - vol_dest
                    if pour_amount > 0:
                        next_state_list = list(current_state)
                        next_state_list[i] -= pour_amount
                        next_state_list[j] = cap_dest
                        next_state = tuple(next_state_list)
                        if next_state not in visited:
                            visited.add(next_state)
                            queue.append((next_state, path + [f"P({name_src}, {name_dest})"]))
                            
    print("No solution found.")

if __name__ == '__main__':
    solve_oil_puzzle()