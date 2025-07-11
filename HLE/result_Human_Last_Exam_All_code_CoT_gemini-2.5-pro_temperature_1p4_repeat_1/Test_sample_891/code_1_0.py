import collections

def solve_oil_division_puzzle():
    """
    Solves the oil pouring puzzle to find the shortest sequence of actions
    to divide 39L of oil into three 13L portions using BFS.
    """
    # Capacities of containers. We use a tuple for easy indexing.
    # Index 0 for X is its initial volume, not a capacity constraint for others.
    capacities = (float('inf'), 8, 17, 21) # (X, A, B, C)
    container_names = ('X', 'A', 'B', 'C')

    # Initial state: (oil_in_x, oil_in_a, oil_in_b, oil_in_c)
    initial_state = (39, 0, 0, 0)
    
    # Target state: Three 13L portions, which must be in X, B, C.
    target_state = (13, 0, 13, 13)
    
    # Queue for BFS. Each item is (state, path_to_state).
    queue = collections.deque([(initial_state, [])])
    
    # A set to keep track of visited states to prevent cycles.
    visited = {initial_state}
    
    while queue:
        current_state, path = queue.popleft()

        # If the target state is reached, we are done.
        if current_state == target_state:
            print("Shortest solution found!")
            print(f"Total steps: {len(path)}")
            # The question asks for the last operation in the sequence.
            last_operation = path[-1]
            print("\nThe final equation is reaching the target state (X:13, A:0, B:13, C:13).")
            # We can find the penultimate state to show the last step clearly.
            # To do so requires re-calculating the entire path.
            # However, the last action itself is sufficient.
            
            # The problem asks to output each number in the final equation.
            # Let's derive the state before the last pour.
            # Last action P(A,B) would mean A was emptied into B.
            # If so, the previous state had A > 0 and B < 13.
            # For example, if the last action is P(A, B) and A had 8L, B must have had 5L.
            # So, 13 (in X) + 8 (in A) + 5 (in B) + 13 (in C) -> P(A, B) -> 13 + 0 + (5+8) + 13
            # which is 13 (in X) + 0 (in A) + 13 (in B) + 13 (in C).
            # The numbers in this final equation are 13, 0, 13, 13.
            final_X, final_A, final_B, final_C = current_state
            print(f"Final amounts: X={final_X}, A={final_A}, B={final_B}, C={final_C}")
            
            print(f"\nThe last operation is: {last_operation}")
            return

        # Generate all valid next states from the current state.
        for from_idx in range(4): # 0:X, 1:A, 2:B, 3:C
            for to_idx in range(1, 4): # Can't pour into X (index 0)
                if from_idx == to_idx:
                    continue

                # Create a mutable list from the state tuple
                next_state_list = list(current_state)
                
                from_amount = next_state_list[from_idx]
                to_amount = next_state_list[to_idx]
                to_capacity = capacities[to_idx]
                
                # Cannot pour from an empty container or into a full one.
                if from_amount == 0 or to_amount == to_capacity:
                    continue

                # Calculate pour amount based on rules (empty source or fill destination)
                amount_to_pour = min(from_amount, to_capacity - to_amount)
                
                if amount_to_pour == 0:
                    continue

                # Perform the pour
                next_state_list[from_idx] -= amount_to_pour
                next_state_list[to_idx] += amount_to_pour

                next_state_tuple = tuple(next_state_list)

                # If we have not visited this state, add it to the queue.
                if next_state_tuple not in visited:
                    visited.add(next_state_tuple)
                    action = f"P({container_names[from_idx]}, {container_names[to_idx]})"
                    new_path = path + [action]
                    queue.append((next_state_tuple, new_path))
    
    print("No solution found.")

# Run the solver.
solve_oil_division_puzzle()