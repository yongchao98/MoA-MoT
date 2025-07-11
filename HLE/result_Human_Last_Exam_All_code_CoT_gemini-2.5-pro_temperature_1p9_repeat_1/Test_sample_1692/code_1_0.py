def solve_cellular_automaton():
    """
    Finds the intermediate step in a 2-step evolution of an elementary cellular automaton.
    
    This function iterates through all 256 possible ECA rules to find the one that
    correctly transforms the initial state to the final state in two steps.
    """
    
    initial_state_str = '01101001'
    final_state_str = '10000111'
    num_cells = len(initial_state_str)
    
    solutions = []

    def apply_rule(state, rule_binary):
        """
        Applies a given ECA rule to a state string and returns the next state.
        Handles periodic (wrap-around) boundary conditions.
        
        The rule_binary is an 8-bit string corresponding to the outputs for
        neighborhoods '111', '110', '101', '100', '011', '010', '001', '000'.
        """
        next_state = []
        # Create a mapping from neighborhood pattern to its index in the rule binary string.
        # This is equivalent to mapping '111' to bit 0, '110' to bit 1, ..., '000' to bit 7.
        rule_map = {format(i, '03b'): rule_binary[7 - i] for i in range(8)}

        for i in range(num_cells):
            # Get neighborhood with periodic boundaries.
            # In Python, state[-1] correctly accesses the last element.
            left = state[i - 1]
            center = state[i]
            # The modulo operator handles the wrap-around for the right neighbor.
            right = state[(i + 1) % num_cells]
            
            neighborhood = left + center + right
            
            # Look up the new cell state from the rule map and append it.
            next_state.append(rule_map[neighborhood])
            
        return "".join(next_state)

    # Iterate through all 256 possible ECA rules.
    for rule_number in range(256):
        # Convert the rule number to its 8-bit binary representation.
        rule_binary_str = format(rule_number, '08b')
        
        # Step 1: Generate the candidate intermediate state from the initial state.
        intermediate_state = apply_rule(initial_state_str, rule_binary_str)
        
        # Step 2: Generate the final state from the candidate intermediate state.
        calculated_final_state = apply_rule(intermediate_state, rule_binary_str)
        
        # Step 3: Compare with the known final state.
        if calculated_final_state == final_state_str:
            solutions.append(intermediate_state)
            # You could break here if you assume there is only one solution,
            # but we'll continue to ensure it's unique.

    # Output the result.
    if len(solutions) == 1:
        # As requested, format the output with spaces between each digit.
        result = " ".join(solutions[0])
        print(f"Initial state:    { ' '.join(initial_state_str) }")
        print(f"Intermediate state: {result}")
        print(f"Final state:        { ' '.join(final_state_str) }")
    elif len(solutions) == 0:
        print("No valid solution found.")
    else:
        print(f"Multiple ({len(solutions)}) solutions found:")
        for sol in solutions:
            print(" ".join(sol))

if __name__ == '__main__':
    solve_cellular_automaton()
<<<01101100>>>