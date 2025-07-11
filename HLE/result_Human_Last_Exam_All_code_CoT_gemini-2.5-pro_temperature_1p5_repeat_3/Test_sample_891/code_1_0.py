import collections

def solve_oil_puzzle():
    """
    Solves the oil pouring puzzle using Breadth-First Search to find the shortest sequence of operations.
    """
    # Capacities of the containers
    capacities = {'X': 39, 'A': 8, 'B': 17, 'C': 21}
    # An ordered list of container names to map indices to names
    container_names = ['X', 'A', 'B', 'C']

    # Initial state: (oil in X, A, B, C)
    initial_state = (39, 0, 0, 0)
    # Target state: three 13L portions
    target_state = (13, 0, 13, 13)

    # Queue for BFS, storing tuples of (current_state, path_to_state)
    queue = collections.deque([(initial_state, [])])
    
    # A set to keep track of visited states to avoid cycles and redundant computations
    visited = {initial_state}

    # All possible pouring moves (from_container, to_container)
    # The rule "You cannot pour back to container X" is implicitly handled
    # by not including moves with 'X' as a destination.
    possible_moves = [
        ('X', 'A'), ('X', 'B'), ('X', 'C'),
        ('A', 'B'), ('A', 'C'),
        ('B', 'A'), ('B', 'C'),
        ('C', 'A'), ('C', 'B'),
    ]

    while queue:
        current_state, path = queue.popleft()

        # If the target state is found, print the solution and exit.
        if current_state == target_state:
            print("Shortest sequence of operations found:")
            
            # Trace the path from start to finish and print each step
            temp_state = list(initial_state)
            print(f"Start: X={temp_state[0]}, A={temp_state[1]}, B={temp_state[2]}, C={temp_state[3]}")
            for i, move_str in enumerate(path):
                source_name, dest_name = move_str[2], move_str[5]
                source_idx = container_names.index(source_name)
                dest_idx = container_names.index(dest_name)
                
                # Calculate how much oil is poured
                amount = min(
                    temp_state[source_idx], 
                    capacities[dest_name] - temp_state[dest_idx]
                )
                
                # Update the state for the printout
                temp_state[source_idx] -= amount
                temp_state[dest_idx] += amount
                
                print(f"{i+1}. P({source_name}, {dest_name}) -> X={temp_state[0]}, A={temp_state[1]}, B={temp_state[2]}, C={temp_state[3]}")

            print("\nFinal State (The 'Final Equation'):")
            final_x, final_a, final_b, final_c = current_state
            print(f"Container X: {final_x}")
            print(f"Container A: {final_a}")
            print(f"Container B: {final_b}")
            print(f"Container C: {final_c}")

            # The question asks for the last operation
            last_operation = path[-1]
            print(f"\nThe last operation to achieve the goal is: {last_operation}")
            return

        # If target not found, generate next possible states
        for source_name, dest_name in possible_moves:
            next_state_list = list(current_state)
            source_idx = container_names.index(source_name)
            dest_idx = container_names.index(dest_name)
            
            # A pour is only possible if the source is not empty and destination is not full
            if next_state_list[source_idx] == 0 or next_state_list[dest_idx] == capacities[dest_name]:
                continue

            # Calculate amount to pour
            amount = min(
                next_state_list[source_idx],
                capacities[dest_name] - next_state_list[dest_idx]
            )

            # Apply the pour to get the next state
            next_state_list[source_idx] -= amount
            next_state_list[dest_idx] += amount
            next_state = tuple(next_state_list)

            # If this state has not been visited, add it to the queue and visited set
            if next_state not in visited:
                visited.add(next_state)
                new_path = path + [f"P({source_name},{dest_name})"]
                queue.append((next_state, new_path))
                
    print("No solution was found.")

if __name__ == '__main__':
    solve_oil_puzzle()