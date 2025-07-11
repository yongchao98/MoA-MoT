from collections import deque

def solve_oil_puzzle():
    """
    Solves the oil pouring puzzle using Breadth-First Search (BFS)
    to find the shortest sequence of operations.
    """
    # Capacities of containers: (X, A, B, C)
    # X is effectively infinite as a source, but we use its initial amount.
    capacities = (39, 8, 17, 21)
    
    # Initial state: (oil_in_X, oil_in_A, oil_in_B, oil_in_C)
    initial_state = (39, 0, 0, 0)
    
    # Target state: three equal portions of 13L.
    # Container A (8L) cannot hold 13L.
    # So, the target is 13L in X, B, and C.
    target_state = (13, 0, 13, 13)

    # Queue for BFS: stores tuples of (state, path)
    queue = deque([(initial_state, [])])
    
    # Set to keep track of visited states to avoid cycles
    visited = {initial_state}

    # Map indices to container names for printing the path
    names = {0: 'X', 1: 'A', 2: 'B', 3: 'C'}

    while queue:
        current_state, path = queue.popleft()

        if current_state == target_state:
            print("Shortest sequence of operations found:")
            # Print the full state transition for clarity
            start_state = list(initial_state)
            print(f"Start: {initial_state}")
            for move in path:
                # This part is for re-calculating the state after each move to print it
                # It's not part of the BFS search itself.
                s_name, d_name = move.replace('P(', '').replace(')', '').split(', ')
                s_idx = [k for k, v in names.items() if v == s_name][0]
                d_idx = [k for k, v in names.items() if v == d_name][0]
                
                s_val = start_state[s_idx]
                d_val = start_state[d_idx]
                d_cap = capacities[d_idx]
                
                amount = min(s_val, d_cap - d_val)
                
                start_state[s_idx] -= amount
                start_state[d_idx] += amount
                
                print(f"-> {move}: {tuple(start_state)}")

            print("\n---")
            print(f"The total number of steps is: {len(path)}")
            print(f"The last operation to achieve the goal is: {path[-1]}")
            return path[-1]

        # Generate all possible next states by pouring between containers
        # i is the source container index, j is the destination container index
        for i in range(4):
            for j in range(1, 4): # Can only pour into A, B, C
                if i == j:
                    continue

                # Create a mutable list from the current state tuple
                next_state_list = list(current_state)
                
                s_val = next_state_list[i]
                d_val = next_state_list[j]
                d_cap = capacities[j]

                # If source is empty or destination is full, skip
                if s_val == 0 or d_val == d_cap:
                    continue

                # Calculate amount to pour
                amount = min(s_val, d_cap - d_val)

                # Update amounts in source and destination
                next_state_list[i] -= amount
                next_state_list[j] += amount
                
                next_state = tuple(next_state_list)

                if next_state not in visited:
                    visited.add(next_state)
                    new_path = path + [f"P({names[i]}, {names[j]})"]
                    queue.append((next_state, new_path))

    print("No solution found.")
    return None

if __name__ == '__main__':
    solve_oil_puzzle()