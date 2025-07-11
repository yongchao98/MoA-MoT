from collections import deque

def solve_pouring_problem():
    """
    Solves the oil pouring puzzle using Breadth-First Search (BFS)
    to find the shortest sequence of operations.
    """
    # Container capacities and names
    capacities = {'A': 8, 'B': 17, 'C': 21}
    names = ['A', 'B', 'C']
    
    # State is represented by a tuple (a, b, c) for the volumes in containers A, B, and C.
    initial_state = (0, 0, 0)
    
    # BFS queue stores tuples of (state, path_to_state)
    queue = deque([(initial_state, [])])
    
    # A set to keep track of visited states to avoid cycles and redundant work.
    visited = {initial_state}

    while queue:
        current_state_tuple, path = queue.popleft()
        a, b, c = current_state_tuple
        
        # Check if the current state meets the goal condition (B=13, C=13)
        if b == 13 and c == 13:
            # Found the solution. Now, print the steps and the final state.
            print("Found the shortest path in", len(path), "steps:")
            
            # Initialize volumes for step-by-step display
            vols = {'X': 39, 'A': 0, 'B': 0, 'C': 0}
            print(f"Start:          (X: {vols['X']}, A: {vols['A']}, B: {vols['B']}, C: {vols['C']})")

            for i, step in enumerate(path):
                source, dest = step.replace("P(", "").replace(")", "").split(", ")
                
                if source == 'X':
                    pour = capacities[dest] - vols[dest]
                    vols[source] -= pour
                    vols[dest] += pour
                else:
                    pour = min(vols[source], capacities[dest] - vols[dest])
                    vols[source] -= pour
                    vols[dest] += pour
                
                print(f"Step {i+1}: {step} -> (X: {vols['X']}, A: {vols['A']}, B: {vols['B']}, C: {vols['C']})")
            
            print("\nFinal distribution:")
            final_x = 39 - a - b - c
            print(f"Portion 1 in Container B: {b}")
            print(f"Portion 2 in Container C: {c}")
            print(f"Portion 3 in (X + A): {final_x} + {a} = {final_x + a}")
            return # Exit after printing the first (and shortest) solution

        # --- Generate next possible states from the current state ---

        # List representation for easier modification
        current_state_list = [a, b, c]
        x_vol = 39 - sum(current_state_list)

        # Action: Pour from source X to fill A, B, or C
        for i in range(3):
            dest_name = names[i]
            # If destination is not full and there's enough oil in X
            if current_state_list[i] < capacities[dest_name] and x_vol >= (capacities[dest_name] - current_state_list[i]):
                next_state_list = list(current_state_list)
                next_state_list[i] = capacities[dest_name]
                next_state = tuple(next_state_list)
                if next_state not in visited:
                    visited.add(next_state)
                    queue.append((next_state, path + [f"P(X, {dest_name})"]))
        
        # Action: Pour between containers A, B, and C
        for i in range(3):  # Source container index
            for j in range(3):  # Destination container index
                if i == j: continue # Cannot pour into itself
                
                pour = min(current_state_list[i], capacities[names[j]] - current_state_list[j])
                if pour > 0:
                    next_state_list = list(current_state_list)
                    next_state_list[i] -= pour
                    next_state_list[j] += pour
                    next_state = tuple(next_state_list)
                    if next_state not in visited:
                        visited.add(next_state)
                        queue.append((next_state, path + [f"P({names[i]}, {names[j]})"]))

solve_pouring_problem()