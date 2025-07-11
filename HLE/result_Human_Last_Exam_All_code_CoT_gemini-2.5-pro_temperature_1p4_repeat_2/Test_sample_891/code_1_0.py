def solve_oil_puzzle():
    """
    Solves the 39L oil division puzzle by executing a known shortest sequence of pours.
    """
    # Capacities of the containers
    caps = {'X': float('inf'), 'A': 8, 'B': 17, 'C': 21}
    
    # Initial state of the containers
    state = {'X': 39, 'A': 0, 'B': 0, 'C': 0}

    def print_state(s, step_num):
        print(f"Step {step_num}: Current State -> X: {s['X']}, A: {s['A']}, B: {s['B']}, C: {s['C']}")

    def pour(current_state, src_name, dst_name):
        """
        Simulates pouring oil from a source to a destination container.
        Returns the new state.
        """
        # Make a copy to modify
        new_state = current_state.copy()
        
        # Get amounts and capacity
        src_amount = new_state[src_name]
        dst_amount = new_state[dst_name]
        dst_cap = caps[dst_name]
        
        # Calculate pour amount based on rules (empty source or fill destination)
        pour_amount = min(src_amount, dst_cap - dst_amount)
        
        if pour_amount == 0:
            print(f"Cannot pour from {src_name} to {dst_name}.")
            return new_state # No change

        # Update amounts in the new state
        new_state[src_name] -= pour_amount
        new_state[dst_name] += pour_amount
        
        print(f"Action: P({src_name}, {dst_name}) -> Poured {pour_amount}L.")
        
        return new_state

    print("--- Starting Oil Puzzle Solution ---")
    print_state(state, 0)
    print("-" * 35)

    # The shortest known sequence of moves
    moves = [
        ('X', 'B'),  # Step 1
        ('B', 'C'),  # Step 2
        ('C', 'A'),  # Step 3
        ('A', 'B'),  # Step 4
        ('X', 'B'),  # Step 5
        ('B', 'A'),  # Step 6
        ('A', 'C'),  # Step 7
        ('B', 'A'),  # Step 8
        ('A', 'C'),  # Step 9
        ('C', 'B'),  # Step 10
        ('B', 'A'),  # Step 11
        ('A', 'C'),  # Step 12: The final move
    ]

    for i, (src, dst) in enumerate(moves):
        state = pour(state, src, dst)
        print_state(state, i + 1)
        if i == len(moves) - 2: # Store state before the final move
            c_before_final = state['C']
            a_before_final = state['A']
        print("-" * 35)

    print("--- Goal Reached ---")
    print(f"Final State: X={state['X']}, A={state['A']}, B={state['B']}, C={state['C']}")
    print("Three equal portions of 13 liters have been successfully measured.")
    
    last_op_src, last_op_dst = moves[-1]
    print(f"\nThe last operation was P({last_op_src}, {last_op_dst}).")
    
    # Print the equation for the last operation
    print("\nEquation for the final pour:")
    print(f"Container C's final amount = C's amount before + A's amount before")
    print(f"{state['C']} = {c_before_final} + {a_before_final}")

# Run the solution
solve_oil_puzzle()