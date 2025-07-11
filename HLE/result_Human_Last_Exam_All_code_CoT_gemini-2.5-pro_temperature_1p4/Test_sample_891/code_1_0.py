import collections

def solve_pouring_problem():
    """
    Solves the oil pouring puzzle using Breadth-First Search (BFS)
    to find the shortest sequence of operations.
    """
    
    # Define capacities for containers A, B, and C.
    # Container X starts with 39L.
    capacities = {
        'X': 39, 'A': 8, 'B': 17, 'C': 21
    }
    
    # The state is represented by a tuple (x, a, b, c)
    initial_state = (39, 0, 0, 0)
    target_state = (13, 0, 13, 13)
    
    # A queue for BFS. Each item is a tuple: (state, path_to_state)
    queue = collections.deque([(initial_state, [])])
    
    # A set to keep track of visited states to prevent cycles.
    visited = {initial_state}
    
    # List of container names for generating the path notation P(Source, Destination)
    container_names = ['X', 'A', 'B', 'C']
    cap_tuple = (capacities['X'], capacities['A'], capacities['B'], capacities['C'])

    while queue:
        current_state, path = queue.popleft()

        # If the target is found, we have the shortest path.
        if current_state == target_state:
            # Print the last operation from the found path.
            if path:
                print(path[-1])
            return

        # Generate all possible next states from the current state.
        # i is the source container index, j is the destination container index.
        for i in range(4): # 0:X, 1:A, 2:B, 3:C
            for j in range(4):
                # Cannot pour from a container to itself.
                if i == j:
                    continue
                
                # Restriction: Cannot pour back into container X.
                if j == 0:
                    continue

                source_val = current_state[i]
                dest_val = current_state[j]
                dest_cap = cap_tuple[j]

                # Cannot pour from an empty container or into a full container.
                if source_val == 0 or dest_val == dest_cap:
                    continue
                
                # Calculate the amount to pour.
                # It's the minimum of what's in the source and the remaining space in the destination.
                amount = min(source_val, dest_cap - dest_val)

                # Create the next state.
                next_state_list = list(current_state)
                next_state_list[i] -= amount
                next_state_list[j] += amount
                next_state = tuple(next_state_list)

                # If we haven't seen this state before, add it to the queue and visited set.
                if next_state not in visited:
                    visited.add(next_state)
                    # Construct the new path.
                    action = f"P({container_names[i]}, {container_names[j]})"
                    new_path = path + [action]
                    queue.append((next_state, new_path))
                    
    print("No solution found.")

# Run the solver to find the last operation.
solve_pouring_problem()
<<<B>>>