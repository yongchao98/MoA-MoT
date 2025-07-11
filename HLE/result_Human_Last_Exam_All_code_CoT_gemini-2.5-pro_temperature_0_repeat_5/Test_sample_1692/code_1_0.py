def solve_cellular_automaton():
    """
    Finds the intermediate step in a 2-step elementary cellular automaton evolution
    by testing all 256 possible rules.
    """
    gen_t0 = "01101001"
    gen_t2 = "10000111"
    n = len(gen_t0)
    
    solution_found = None

    # Iterate through all 256 possible rules
    for rule_num in range(256):
        # Get the 8-bit binary representation of the rule
        # This string maps neighborhoods '111', '110', ..., '000' to their output
        rule_binary = format(rule_num, '08b')
        rule_map = {
            '111': rule_binary[0], '110': rule_binary[1], '101': rule_binary[2],
            '100': rule_binary[3], '011': rule_binary[4], '010': rule_binary[5],
            '001': rule_binary[6], '000': rule_binary[7]
        }

        # --- Step 1: Evolve from gen_t0 to gen_t1 ---
        gen_t1 = ""
        for i in range(n):
            # Get neighborhood with periodic (wrap-around) boundaries
            left = gen_t0[i - 1]
            center = gen_t0[i]
            right = gen_t0[(i + 1) % n]
            neighborhood = left + center + right
            gen_t1 += rule_map[neighborhood]

        # --- Step 2: Evolve from gen_t1 to a calculated gen_t2 ---
        calculated_gen_t2 = ""
        for i in range(n):
            left = gen_t1[i - 1]
            center = gen_t1[i]
            right = gen_t1[(i + 1) % n]
            neighborhood = left + center + right
            calculated_gen_t2 += rule_map[neighborhood]

        # --- Step 3: Check if the calculated final state matches the given one ---
        if calculated_gen_t2 == gen_t2:
            solution_found = gen_t1
            break
    
    # The problem guarantees a unique solution is found
    if solution_found:
        print(solution_found)

solve_cellular_automaton()