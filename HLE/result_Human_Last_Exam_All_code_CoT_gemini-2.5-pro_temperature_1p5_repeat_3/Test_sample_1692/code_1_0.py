def solve_cellular_automaton():
    """
    Finds the intermediate step in a 1D elementary cellular automaton sequence
    by testing all 256 possible rules.
    """

    gen1 = "01101001"
    gen3 = "10000111"
    width = len(gen1)
    
    # Helper function to apply a rule to a state
    def apply_rule(state, rule_num):
        # Convert the rule number to its 8-bit binary representation.
        # This binary string defines the output for the 8 possible neighborhoods.
        # e.g., Rule 30 is 00011110 in binary. The '1' at index 3 means '100' -> 1.
        rule_bin = format(rule_num, '08b')
        rule_map = {
            (1, 1, 1): int(rule_bin[0]),
            (1, 1, 0): int(rule_bin[1]),
            (1, 0, 1): int(rule_bin[2]),
            (1, 0, 0): int(rule_bin[3]),
            (0, 1, 1): int(rule_bin[4]),
            (0, 1, 0): int(rule_bin[5]),
            (0, 0, 1): int(rule_bin[6]),
            (0, 0, 0): int(rule_bin[7]),
        }

        # Convert state string to a list of ints
        state_int = [int(c) for c in state]
        next_state_int = []

        for i in range(width):
            # Get neighbors using periodic boundary conditions.
            # Python's negative index `state_int[-1]` handles the left boundary for i=0.
            left = state_int[i - 1]
            center = state_int[i]
            # Modulo operator `%` handles the right boundary for i=width-1.
            right = state_int[(i + 1) % width]
            
            neighborhood = (left, center, right)
            next_state_int.append(rule_map[neighborhood])

        # Convert list of ints back to a string
        return "".join(map(str, next_state_int))

    # Iterate through all 256 possible rules
    for rule in range(256):
        # Generate the intermediate state from generation 1
        gen2_candidate = apply_rule(gen1, rule)
        # Generate the next state from the intermediate state
        gen3_candidate = apply_rule(gen2_candidate, rule)

        # Check if the generated state matches the known generation 3
        if gen3_candidate == gen3:
            # We found the unique solution. Print each digit of the intermediate state.
            # The puzzle is:
            # 0 1 1 0 1 0 0 1
            # ? ? ? ? ? ? ? ?  <-- This is what we are printing
            # 1 0 0 0 0 1 1 1
            print(' '.join(list(gen2_candidate)))
            return gen2_candidate

# Run the solver
solution = solve_cellular_automaton()
print(f"<<<{solution}>>>")