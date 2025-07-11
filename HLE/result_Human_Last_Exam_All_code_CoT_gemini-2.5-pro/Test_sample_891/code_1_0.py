def solve_oil_puzzle():
    """
    This function simulates the steps to solve the oil division puzzle
    and identifies the final operation.
    """

    # Initial state of the containers in liters
    state = {'X': 39, 'A': 0, 'B': 0, 'C': 0}
    # Capacities of the containers
    caps = {'A': 8, 'B': 17, 'C': 21}

    def pour(current_state, source_key, dest_key):
        """Simulates a single pour operation."""
        new_state = current_state.copy()
        
        # Pouring from X means filling the destination
        if source_key == 'X':
            space_in_dest = caps[dest_key] - new_state[dest_key]
            amount = min(new_state[source_key], space_in_dest)
        # Pouring between A, B, C means emptying the source or filling the destination
        else:
            space_in_dest = caps[dest_key] - new_state[dest_key]
            amount = min(new_state[source_key], space_in_dest)
        
        new_state[source_key] -= amount
        new_state[dest_key] += amount
        return new_state

    # The 7-step solution sequence
    sequence = [
        ('X', 'B'), # 1
        ('B', 'C'), # 2
        ('X', 'B'), # 3
        ('B', 'A'), # 4
        ('A', 'C'), # 5
        ('C', 'B'), # 6
        ('B', 'A'), # 7
    ]

    print("Simulating the pouring sequence:")
    print(f"Start: X={state['X']}, A={state['A']}, B={state['B']}, C={state['C']}")

    for i, (s, d) in enumerate(sequence):
        state = pour(state, s, d)
        print(f"Step {i+1}: P({s}, {d}) -> X={state['X']}, A={state['A']}, B={state['B']}, C={state['C']}")

    print("\n--- Final State Analysis ---")
    print(f"Portion 1: Container B has {state['B']} liters.")
    print(f"Portion 2: Container C has {state['C']} liters.")
    print(f"Portion 3: Containers X and A together have {state['X']} + {state['A']} = {state['X'] + state['A']} liters.")

    last_op_s, last_op_d = sequence[-1]
    last_operation_string = f"P({last_op_s}, {last_op_d})"
    
    print(f"\nThe goal is achieved. The last operation in the sequence is: {last_operation_string}")

solve_oil_puzzle()