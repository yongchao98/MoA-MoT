from collections import deque

def solve_pouring_puzzle():
    """
    Solves the 39L oil pouring puzzle using Breadth-First Search (BFS)
    to find the shortest sequence of operations.
    """
    # Capacities of the containers A, B, C
    capacities = {'A': 8, 'B': 17, 'C': 21}
    container_names = ['A', 'B', 'C']
    
    # Initial state: (amount_A, amount_B, amount_C)
    initial_state = (0, 0, 0)
    
    # Target states: Two containers have 13L, one is empty.
    # This implies the main container X also holds 13L.
    target_states = {
        (13, 13, 0),
        (13, 0, 13),
        (0, 13, 13)
    }
    
    # The queue for BFS will store tuples of (state, path)
    queue = deque([(initial_state, [])])
    
    # A set to keep track of visited states to avoid cycles and redundant work
    visited = {initial_state}
    
    while queue:
        current_state, path = queue.popleft()
        
        # Check if we have reached one of the target states
        if current_state in target_states:
            final_a, final_b, final_c = current_state
            final_x = 39 - (final_a + final_b + final_c)
            
            print(f"Found a shortest solution in {len(path)} steps.")
            print("\n--- Sequence of Operations ---")
            
            # To display the state at each step, we re-calculate the trace
            s = {'A': 0, 'B': 0, 'C': 0}
            x = 39
            print(f"Start: X={x:2d}, A={s['A']:2d}, B={s['B']:2d}, C={s['C']:2d}")
            for move in path:
                source, dest = move.replace("P(", "").replace(")", "").split(", ")
                
                s_amount = x if source == 'X' else s[source]
                d_amount = s[dest]
                d_cap = capacities[dest]
                
                pour = min(s_amount, d_cap - d_amount)
                
                if source == 'X':
                    x -= pour
                else:
                    s[source] -= pour
                s[dest] += pour
                print(f"{move:>7s}: X={x:2d}, A={s['A']:2d}, B={s['B']:2d}, C={s['C']:2d}")

            print("\n--- Final Equation ---")
            print("The final distribution of oil is:")
            print(f"Container X = {final_x} L")
            print(f"Container A = {final_a} L")
            print(f"Container B = {final_b} L")
            print(f"Container C = {final_c} L")
            
            print("\n--- Answer ---")
            print(f"The last operation to achieve the goal is: {path[-1]}")
            return
            
        # Generate next possible states from the current state
        a, b, c = current_state
        current_amounts = {'A': a, 'B': b, 'C': c}
        x = 39 - (a + b + c)
        
        all_sources = {'X': x, 'A': a, 'B': b, 'C': c}
        
        # Iterate through all possible pouring actions (Source -> Destination)
        for s_name, s_amount in all_sources.items():
            if s_amount == 0:
                continue
                
            for d_name in container_names:
                if s_name == d_name:
                    continue
                
                d_amount = current_amounts[d_name]
                d_capacity = capacities[d_name]
                
                if d_amount < d_capacity:
                    pour_amount = min(s_amount, d_capacity - d_amount)
                    
                    next_amounts = list(current_state)
                    if s_name != 'X':
                        s_index = container_names.index(s_name)
                        next_amounts[s_index] -= pour_amount
                    
                    d_index = container_names.index(d_name)
                    next_amounts[d_index] += pour_amount
                    next_state = tuple(next_amounts)
                    
                    if next_state not in visited:
                        visited.add(next_state)
                        new_path = path + [f"P({s_name}, {d_name})"]
                        queue.append((next_state, new_path))
                        
    print("No solution found.")

# Execute the solver
solve_pouring_puzzle()