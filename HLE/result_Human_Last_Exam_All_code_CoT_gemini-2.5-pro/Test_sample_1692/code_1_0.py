def solve_cellular_automaton():
    """
    Finds the unique intermediate step in a 1D elementary cellular automaton
    by testing all 256 possible rules.
    """
    gen_t0_str = "01101001"
    gen_t2_str = "10000111"

    gen_t0 = [int(c) for c in gen_t0_str]
    length = len(gen_t0)

    solutions = []

    def apply_rule(rule_num, generation):
        """
        Applies a given rule to a generation to produce the next generation.
        Assumes periodic boundary conditions.
        """
        # Convert rule number to 8-bit binary string (e.g., 30 -> "00011110")
        rule_bits = format(rule_num, '08b')
        
        # Map the 8 neighborhood patterns (from '111' down to '000') to the rule's output
        rule_map = {
            (1, 1, 1): int(rule_bits[0]), (1, 1, 0): int(rule_bits[1]),
            (1, 0, 1): int(rule_bits[2]), (1, 0, 0): int(rule_bits[3]),
            (0, 1, 1): int(rule_bits[4]), (0, 1, 0): int(rule_bits[5]),
            (0, 0, 1): int(rule_bits[6]), (0, 0, 0): int(rule_bits[7])
        }

        next_gen = [0] * length
        for i in range(length):
            # Get neighbors with periodic (wrap-around) boundaries
            left = generation[i - 1]
            center = generation[i]
            right = generation[(i + 1) % length]
            
            neighborhood = (left, center, right)
            next_gen[i] = rule_map[neighborhood]
            
        return next_gen

    # Iterate through all 256 possible rules
    for rule in range(256):
        # Calculate the potential intermediate generation from the first generation
        gen_t1_candidate = apply_rule(rule, gen_t0)

        # Calculate the third generation from the candidate intermediate generation
        gen_t2_calculated = apply_rule(rule, gen_t1_candidate)

        # Convert the result back to a string for comparison
        gen_t2_calculated_str = "".join(map(str, gen_t2_calculated))

        # Check if the calculated third generation matches the given one
        if gen_t2_calculated_str == gen_t2_str:
            solution_str = "".join(map(str, gen_t1_candidate))
            if solution_str not in solutions:
                solutions.append(solution_str)

    # As per the problem, there is only one valid solution. Print it.
    if len(solutions) == 1:
        # The prompt asks to "output each number in the final equation!"
        # This is interpreted as printing the digits of the solution separated by spaces.
        print(" ".join(solutions[0]))
    elif len(solutions) > 1:
        print("Multiple solutions found:", solutions)
    else:
        print("No solution found.")

solve_cellular_automaton()