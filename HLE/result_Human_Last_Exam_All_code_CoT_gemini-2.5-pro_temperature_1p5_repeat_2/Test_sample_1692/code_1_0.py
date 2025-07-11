def find_intermediate_step():
    """
    Solves for the intermediate state of an elementary cellular automaton
    by iterating through all 256 possible rules.
    """
    gen0 = "01101001"
    gen2_target = "10000111"

    def apply_rule(state, rule_num):
        """
        Applies a given cellular automaton rule to a state, assuming periodic boundaries.
        
        Args:
            state (str): The binary string representing the current state.
            rule_num (int): The rule number (0-255).
        
        Returns:
            str: The binary string representing the next state.
        """
        # Convert the rule number to an 8-bit binary string representation.
        # The bits correspond to the output for neighborhoods '111' through '000'.
        rule_bits = format(rule_num, '08b')
        
        # Create a mapping from neighborhood pattern to output bit
        patterns = {
            '111': rule_bits[0], '110': rule_bits[1], '101': rule_bits[2],
            '100': rule_bits[3], '011': rule_bits[4], '010': rule_bits[5],
            '001': rule_bits[6], '000': rule_bits[7]
        }
        
        n = len(state)
        next_state = []
        for i in range(n):
            # Get the 3-cell neighborhood using periodic boundary conditions.
            left = state[(i - 1 + n) % n]
            center = state[i]
            right = state[(i + 1) % n]
            neighborhood = left + center + right
            
            # Determine the new state for the cell and append it.
            next_state.append(patterns[neighborhood])
            
        return "".join(next_state)

    # Iterate through all 256 possible elementary rules.
    for rule in range(256):
        # First step: gen0 -> gen1
        gen1_candidate = apply_rule(gen0, rule)
        
        # Second step: gen1 -> gen2
        gen2_candidate = apply_rule(gen1_candidate, rule)
        
        # Check if the generated final state matches the target.
        if gen2_candidate == gen2_target:
            # If it matches, we have found the solution.
            print(gen1_candidate)
            return

# Execute the function to find and print the solution.
find_intermediate_step()