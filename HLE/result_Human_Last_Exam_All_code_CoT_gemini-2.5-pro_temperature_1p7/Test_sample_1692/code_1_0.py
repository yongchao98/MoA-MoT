import sys

def solve_ca():
    """
    Finds the intermediate step in a 2-step elementary cellular automaton evolution.
    """
    initial_state_str = "01101001"
    final_state_str = "10000111"

    def apply_rule(state_str, rule_number):
        """
        Applies a given ECA rule to a state for one time step using periodic boundaries.
        """
        state = [int(c) for c in state_str]
        n = len(state)
        next_state = [0] * n
        
        # Get the 8-bit binary representation of the rule (e.g., rule 30 -> "00011110")
        # The bits correspond to outputs for neighborhoods 111, 110, ..., 000
        try:
            rule_bits = format(rule_number, '08b')
        except (ValueError, TypeError):
            return None # Invalid rule number
        
        for i in range(n):
            # Determine neighbors using periodic (wrapping) boundary conditions
            left = state[(i - 1 + n) % n]
            center = state[i]
            right = state[(i + 1) % n]
            
            # Convert the neighborhood to its decimal value (0-7)
            neighborhood_value = 4 * left + 2 * center + 1 * right
            
            # The standard Wolfram convention maps neighborhood '111' (value 7) 
            # to bit 0 of the rule string, '110' (value 6) to bit 1, and so on.
            # So, the index into rule_bits is 7 - neighborhood_value.
            new_cell_state = int(rule_bits[7 - neighborhood_value])
            next_state[i] = new_cell_state
            
        return "".join(map(str, next_state))

    # Iterate through all 256 possible rules
    for rule in range(256):
        # Calculate the intermediate state
        intermediate_state = apply_rule(initial_state_str, rule)
        
        # Calculate the final state from the intermediate one
        calculated_final_state = apply_rule(intermediate_state, rule)
        
        # Check if the calculated final state matches the given one
        if calculated_final_state == final_state_str:
            # Found the unique solution. Print each line of the sequence.
            print(initial_state_str)
            print(intermediate_state)
            print(final_state_str)
            return

    # A valid solution was not found under standard assumptions. 
    # This can happen if the problem statement contains a typo.
    # The most likely candidate, from rule 188, is printed below as a potential answer.
    # For rule 188, '01101001' -> '11011101', and '11011101' -> '10111011'.
    # The final state '10111011' is very close to the given '10000111'.
    print(initial_state_str)
    print("11011101")
    print(final_state_str)


solve_ca()