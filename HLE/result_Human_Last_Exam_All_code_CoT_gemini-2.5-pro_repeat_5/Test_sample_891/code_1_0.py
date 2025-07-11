from collections import deque

def solve_pouring_problem():
    """
    Solves the oil pouring puzzle using Breadth-First Search (BFS)
    to find the shortest sequence of operations.
    """
    # Capacities of containers (A, B, C)
    capacities = (8, 17, 21)
    # Target state for containers (A, B, C)
    target_state = (0, 13, 13)

    # The queue for BFS will store tuples of (state, path)
    # state is a tuple (a, b, c) of oil amounts
    # path is a list of tuples, each containing the action string and the resulting state
    queue = deque([(((0, 0, 0), [("Initial State", (0, 0, 0))]))])

    # A set to keep track of visited states to avoid cycles
    visited = {(0, 0, 0)}

    while queue:
        (a, b, c), path = queue.popleft()

        # If the target is reached, print the solution and return
        if (a, b, c) == target_state:
            print("Solution found!")
            print("Step-by-step sequence:")
            
            initial_X = 39
            for i, (action, (ca, cb, cc)) in enumerate(path):
                oil_in_X = initial_X - (ca + cb + cc)
                print(f"{i}. {action:<15} -> State (A,B,C): ({ca}, {cb}, {cc}), X remaining: {oil_in_X}")

            final_equation_parts = []
            for i, (action, (ca, cb, cc)) in enumerate(path):
                if i > 0:
                    final_equation_parts.append(action)

            # To satisfy the output format requirement of showing the final equation.
            print("\nFinal sequence of operations:")
            # The problem asks for the last operation in the equation.
            # We will print all operations in the sequence for clarity.
            final_action = path[-1][0]
            for part in final_equation_parts:
                print(part)

            print(f"\nThe last operation is: {final_action}")
            return

        # --- Generate next possible states ---
        
        # Current amount of oil in containers and in X
        current_sum = a + b + c
        oil_in_x = 39 - current_sum
        
        # Define sources and destinations for pouring
        # Format: (source_name, source_val, dest_name, dest_val, dest_cap)
        # Note: 'X' is a source, but not a destination
        
        # Pouring from X
        sources = {'A': a, 'B': b, 'C': c}
        for dest_name, dest_val in sources.items():
            dest_idx = ord(dest_name) - ord('A')
            dest_cap = capacities[dest_idx]
            
            if dest_val < dest_cap:
                amount_needed = dest_cap - dest_val
                if oil_in_x >= amount_needed:
                    next_state_list = list((a, b, c))
                    next_state_list[dest_idx] = dest_cap
                    next_state = tuple(next_state_list)
                    if next_state not in visited:
                        visited.add(next_state)
                        action_str = f"P(X, {dest_name})"
                        queue.append((next_state, path + [(action_str, next_state)]))
        
        # Pouring between A, B, C
        containers = [('A', a, capacities[0]), ('B', b, capacities[1]), ('C', c, capacities[2])]
        for i in range(3):
            for j in range(3):
                if i == j:
                    continue
                
                s_name, s_val, s_cap = containers[i]
                d_name, d_val, d_cap = containers[j]
                
                pour_amount = min(s_val, d_cap - d_val)
                
                if pour_amount > 0:
                    next_state_list = list((a, b, c))
                    next_state_list[i] -= pour_amount
                    next_state_list[j] += pour_amount
                    next_state = tuple(next_state_list)
                    
                    if next_state not in visited:
                        visited.add(next_state)
                        action_str = f"P({s_name}, {d_name})"
                        queue.append((next_state, path + [(action_str, next_state)]))

solve_pouring_problem()
<<<A>>>