import collections

def solve_pouring_problem():
    """
    Solves the oil pouring problem using Breadth-First Search (BFS)
    to find the shortest sequence of operations.
    """
    # Capacities of containers A, B, C. X is considered to have infinite capacity.
    capacities = {'X': float('inf'), 'A': 8, 'B': 17, 'C': 21}
    # For state representation using tuples (X, A, B, C)
    capacities_tuple = (float('inf'), 8, 17, 21)
    container_names = ('X', 'A', 'B', 'C')

    # Initial and target states
    initial_state = (39, 0, 0, 0)
    target_state = (13, 0, 13, 13)

    # Queue for BFS. Each item is a tuple: (state, path_list)
    # The path_list stores tuples of (action_string, next_state)
    queue = collections.deque([(initial_state, [])])
    # A set to store visited states to avoid cycles and redundant work
    visited = {initial_state}

    while queue:
        current_state, path = queue.popleft()

        if current_state == target_state:
            print("Shortest sequence found:")
            start_state_str = f"Start: (X={initial_state[0]}, A={initial_state[1]}, B={initial_state[2]}, C={initial_state[3]})"
            print(start_state_str)

            for action, state in path:
                state_str = f"(X={state[0]}, A={state[1]}, B={state[2]}, C={state[3]})"
                print(f"{action}: {state_str}")

            # The final equation is the target state itself.
            final_equation = f"Final State: X={target_state[0]}, A={target_state[1]}, B={target_state[2]}, C={target_state[3]}"
            # The sum is 13+13+13=39, which is the total oil.
            print(f"\nThis results in three portions of 13 liters in containers X, B, and C.")
            
            # Extract the last operation to answer the multiple-choice question.
            last_operation_text = path[-1][0]
            answer_map = {
                "P(A, B)": "A", "P(A, C)": "B", "P(B, C)": "C",
                "P(C, B)": "D", "P(C, A)": "E", "P(X, A)": "F",
                "P(X, B)": "G", "P(B, A)": "H"
            }
            final_answer = answer_map.get(last_operation_text, "Unknown")
            print(f"<<<{final_answer}>>>")
            return

        # Generate all possible next states from the current state
        for i in range(len(container_names)):  # Source container index
            for j in range(len(container_names)):  # Destination container index
                if i == j:
                    continue

                # Restriction: Cannot pour back to container X (index 0)
                if j == 0:
                    continue

                source_oil = current_state[i]
                dest_oil = current_state[j]
                dest_capacity = capacities_tuple[j]

                # Can't pour from an empty container or into a full one
                if source_oil == 0 or dest_oil == dest_capacity:
                    continue

                # Calculate amount to pour
                pour_amount = min(source_oil, dest_capacity - dest_oil)

                if pour_amount > 0:
                    next_state_list = list(current_state)
                    next_state_list[i] -= pour_amount
                    next_state_list[j] += pour_amount
                    next_state = tuple(next_state_list)

                    if next_state not in visited:
                        visited.add(next_state)
                        action_str = f"P({container_names[i]}, {container_names[j]})"
                        new_path = path + [(action_str, next_state)]
                        queue.append((next_state, new_path))
                        
    print("No solution found.")

solve_pouring_problem()