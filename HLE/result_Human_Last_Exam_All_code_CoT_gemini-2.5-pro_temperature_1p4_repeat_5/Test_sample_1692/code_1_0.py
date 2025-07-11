def solve_cellular_automaton():
    """
    Finds the intermediate step in a 2-generation elementary cellular automaton evolution.
    
    This function iterates through all 256 possible ECA rules, applies each one twice
    to the initial state, and checks if the result matches the final state.
    It assumes periodic boundary conditions.
    """
    initial_str = "01101001"
    final_str = "10000111"

    initial_state = [int(c) for c in initial_str]
    final_state = [int(c) for c in final_str]
    n = len(initial_state)

    # This helper function applies a given rule to a state.
    def apply_rule(current_state, rule_num):
        # The rule number (0-255) defines the output for the 8 possible neighborhoods.
        # Its 8-bit binary representation (e.g., rule 30 is '00011110') maps to the
        # neighborhoods from '111' down to '000'.
        rule_bin = format(rule_num, '08b')
        next_state = []
        
        for i in range(n):
            # Get the 3-cell neighborhood, wrapping around for the edges (periodic boundary).
            left_neighbor = current_state[(i - 1 + n) % n]
            current_cell = current_state[i]
            right_neighbor = current_state[(i + 1) % n]
            
            neighborhood_str = f"{left_neighbor}{current_cell}{right_neighbor}"
            
            # Convert the binary neighborhood string to its integer value (e.g., '111' -> 7).
            neighborhood_idx = int(neighborhood_str, 2)
            
            # The rule's output bit is found by this mapping. '111' (idx 7) maps to rule_bin[0],
            # '110' (idx 6) to rule_bin[1], ..., '000' (idx 0) to rule_bin[7].
            new_cell_val = int(rule_bin[7 - neighborhood_idx])
            next_state.append(new_cell_val)
            
        return next_state

    # Iterate through all 256 possible rules.
    for rule in range(256):
        # Step 1: Generate the intermediate state.
        intermediate_state = apply_rule(initial_state, rule)
        
        # Step 2: Generate the final state from the intermediate one.
        calculated_final = apply_rule(intermediate_state, rule)
        
        # Step 3: Check if the calculated final state matches the target.
        if calculated_final == final_state:
            intermediate_str = "".join(map(str, intermediate_state))
            # If a match is found, print the complete sequence and stop.
            print(initial_str)
            print(intermediate_str)
            print(final_str)
            return

# Execute the solver function.
solve_cellular_automaton()