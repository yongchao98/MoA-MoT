def solve_automaton_puzzle():
    """
    Finds the intermediate step in a 1D elementary cellular automaton evolution
    by testing all 256 possible rules.
    """
    initial_state = "01101001"
    target_final_state = "10000111"

    def apply_rule(rule_num, gen_str):
        """
        Applies a given ECA rule to a generation string to produce the next generation.
        """
        # Per Wolfram convention, rule bits correspond to neighborhood integer values.
        # We reverse the binary string so rule_bits[0] corresponds to '000', etc.
        rule_bits = [int(c) for c in reversed(format(rule_num, '08b'))]
        
        gen_int = [int(c) for c in gen_str]
        n = len(gen_int)
        next_gen_int = [0] * n

        for i in range(n):
            # Get neighborhood values with periodic (wraparound) boundary conditions.
            left_val = gen_int[(i - 1 + n) % n]
            center_val = gen_int[i]
            right_val = gen_int[(i + 1) % n]
            
            # The neighborhood's integer value (0-7) is the index for the rule bit.
            neighborhood_index = 4 * left_val + 2 * center_val + 1 * right_val
            next_gen_int[i] = rule_bits[neighborhood_index]
            
        return "".join(map(str, next_gen_int))

    # Iterate through all 256 possible rules.
    solutions_found = []
    for rule in range(256):
        intermediate_state = apply_rule(rule, initial_state)
        final_state = apply_rule(rule, intermediate_state)
        
        if final_state == target_final_state:
            # A valid solution is found.
            solutions_found.append(intermediate_state)
            
    # As per the prompt, a unique solution is expected.
    if len(solutions_found) == 1:
        the_solution = solutions_found[0]
        # Print the full evolution sequence as requested.
        print("The only valid solution follows the equation:")
        print(f"{initial_state} -> {the_solution} -> {target_final_state}")
    elif len(solutions_found) == 0:
        print("No valid solution was found.")
    else:
        print(f"Found {len(solutions_found)} solutions, expected only one:")
        for sol in solutions_found:
            print(sol)

# Execute the solver function.
solve_automaton_puzzle()