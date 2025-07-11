from collections import deque

def solve_pouring_problem():
    """
    Solves the oil pouring puzzle using Breadth-First Search (BFS)
    to find the shortest sequence of operations.
    """
    # Capacities of containers (X, A, B, C)
    # None for X as its initial volume is the only limit
    capacities = (None, 8, 17, 21)
    
    # Initial and target states
    initial_state = (39, 0, 0, 0)
    # The goal is three 13L portions. Given capacities, this must be
    # 13L in X, 13L in B, and 13L in C.
    goal_state = (13, 0, 13, 13)

    # Queue for BFS: stores tuples of (state, path_list)
    queue = deque([(initial_state, [])])
    # Set to keep track of visited states to avoid cycles
    visited = {initial_state}
    
    container_names = ('X', 'A', 'B', 'C')

    # Start BFS
    while queue:
        current_state, path = queue.popleft()

        # Check if we have reached the goal
        if current_state == goal_state:
            print("Goal reached!")
            print(f"Initial State: (39, 0, 0, 0)")
            
            temp_state = list(initial_state)
            # Print each step and the resulting state
            for i, step in enumerate(path):
                source_name, dest_name = step.replace('P(', '').replace(')', '').split(', ')
                
                # To show the equation, we recalculate the poured amount
                s_idx = container_names.index(source_name)
                d_idx = container_names.index(dest_name)
                
                # Amount to pour
                oil_in_source = temp_state[s_idx]
                if s_idx == 0: # Source is X
                    space_in_dest = capacities[d_idx] - temp_state[d_idx]
                    poured_amount = min(oil_in_source, space_in_dest)
                else: # Source is A, B, or C
                    space_in_dest = capacities[d_idx] - temp_state[d_idx]
                    poured_amount = min(oil_in_source, space_in_dest)

                # Store the equation numbers
                s_before = temp_state[s_idx]
                d_before = temp_state[d_idx]

                # Update state
                temp_state[s_idx] -= poured_amount
                temp_state[d_idx] += poured_amount

                # Print the formatted step
                print(f"Step {i+1}: {step}")
                print(f"  Container {source_name}: {s_before} L -> {temp_state[s_idx]} L")
                print(f"  Container {dest_name}: {d_before} L -> {temp_state[d_idx]} L")
                print(f"  Resulting State: {tuple(temp_state)}")

            print("\nFinal solution path found:")
            print(path)
            
            last_operation = path[-1]
            print(f"\nThe last operation is: {last_operation}")
            return last_operation

        # Generate next possible states
        current_vols = list(current_state)

        # Iterate through all possible pours (Source -> Destination)
        for s_idx, source_name in enumerate(container_names):
            if current_vols[s_idx] == 0:
                continue # Cannot pour from an empty container

            for d_idx, dest_name in enumerate(container_names[1:], 1): # Can't pour to X
                if s_idx == d_idx:
                    continue # Cannot pour into the same container
                
                # Destination is full
                if current_vols[d_idx] == capacities[d_idx]:
                    continue

                next_vols = list(current_vols)
                
                # Amount to pour is the minimum of what's in the source
                # and the available space in the destination.
                if s_idx == 0: # Source is X, which doesn't have a fixed capacity other than current oil
                   space_in_dest = capacities[d_idx] - next_vols[d_idx]
                   poured_amount = min(next_vols[s_idx], space_in_dest)
                else:
                   space_in_dest = capacities[d_idx] - next_vols[d_idx]
                   poured_amount = min(next_vols[s_idx], space_in_dest)
                
                if poured_amount == 0:
                    continue

                # Apply the pour
                next_vols[s_idx] -= poured_amount
                next_vols[d_idx] += poured_amount
                
                next_state = tuple(next_vols)

                if next_state not in visited:
                    visited.add(next_state)
                    new_path = path + [f"P({source_name}, {dest_name})"]
                    queue.append((next_state, new_path))
                    
    return None # Should not be reached if a solution exists

if __name__ == '__main__':
    solution = solve_pouring_problem()
    if not solution:
        print("No solution found.")

<<<A>>>