def solve_cellular_automaton():
    """
    Finds the intermediate step in a 2-step elementary cellular automaton evolution.
    """
    initial_state_str = "01101001"
    final_state_str = "10000111"

    initial_state = [int(c) for c in initial_state_str]
    final_state = [int(c) for c in final_state_str]
    n = len(initial_state)

    def apply_rule(state, rule_map):
        """
        Applies a given rule to a state to get the next state, assuming periodic boundaries.
        """
        next_state = [0] * n
        for i in range(n):
            # Get neighbors with periodic (wrapping) boundary conditions
            left = state[(i - 1 + n) % n]
            center = state[i]
            right = state[(i + 1) % n]
            
            # The 8 neighborhoods from '111' to '000' correspond to indices 0 to 7 in the rule map.
            # Convert the binary neighborhood triplet to an integer value.
            neighborhood_val = left * 4 + center * 2 + right * 1
            
            # The index in the rule map is 7 minus this value.
            rule_index = 7 - neighborhood_val
            
            next_state[i] = rule_map[rule_index]
        return next_state

    # Iterate through all 256 possible elementary rules
    for rule_number in range(256):
        # The rule map is the 8-bit binary representation of the rule number
        rule_map = [int(bit) for bit in format(rule_number, '08b')]
        
        # First step: from initial to intermediate
        intermediate_state = apply_rule(initial_state, rule_map)
        
        # Second step: from intermediate to a calculated final
        calculated_final_state = apply_rule(intermediate_state, rule_map)
        
        # Check if the calculated final state matches the given final state
        if calculated_final_state == final_state:
            # Found the solution. Print the intermediate state.
            solution_str = "".join(map(str, intermediate_state))
            print(solution_str)
            return

solve_cellular_automaton()
<<<00001111>>>