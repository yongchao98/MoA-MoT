def solve_cellular_automaton():
    """
    Finds the intermediate step in a 1D elementary cellular automaton sequence
    by testing all 256 possible rules.
    """
    initial_state = "01101001"
    final_state = "10000111"

    def get_next_generation(state, rule_number):
        """
        Calculates the next generation of a 1D cellular automaton state
        for a given rule, assuming periodic boundary conditions.
        """
        # Convert rule number to its 8-bit binary representation.
        # This string maps outcomes to neighborhoods from '111' down to '000'.
        rule_bits = format(rule_number, '08b')
        
        # Create a mapping from neighborhood patterns to outcomes
        neighborhoods = {
            "111": rule_bits[0], "110": rule_bits[1], "101": rule_bits[2],
            "100": rule_bits[3], "011": rule_bits[4], "010": rule_bits[5],
            "001": rule_bits[6], "000": rule_bits[7]
        }
        
        n = len(state)
        next_state = []
        for i in range(n):
            # Get the 3-cell neighborhood, wrapping around at the edges
            left = state[(i - 1 + n) % n]
            center = state[i]
            right = state[(i + 1) % n]
            
            pattern = f"{left}{center}{right}"
            
            # Append the corresponding outcome to the next state
            next_state.append(neighborhoods[pattern])
            
        return "".join(next_state)

    # Iterate through all 256 possible ECA rules
    for rule in range(256):
        # Step 1: initial_state -> intermediate_state
        intermediate_state = get_next_generation(initial_state, rule)
        
        # Step 2: intermediate_state -> calculated_final_state
        calculated_final_state = get_next_generation(intermediate_state, rule)
        
        # Check if we found the sequence that matches the given final state
        if calculated_final_state == final_state:
            # Print the solution for the intermediate step
            print(intermediate_state)
            return

# Execute the solver
solve_cellular_automaton()