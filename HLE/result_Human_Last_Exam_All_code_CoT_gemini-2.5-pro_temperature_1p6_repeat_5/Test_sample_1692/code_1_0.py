def solve_automaton():
    """
    Finds the intermediate step in a 1D elementary cellular automaton sequence
    by testing all 256 possible rules.
    """
    gen0_str = "01101001"
    gen2_str = "10000111"

    def evolve(state_str, rule_num):
        """
        Calculates the next generation of a cellular automaton.
        Uses periodic boundary conditions.
        """
        # Convert rule number to its 8-bit binary representation.
        # This string directly maps neighborhoods (from '111' down to '000') to outputs.
        rule_bin = format(rule_num, '08b')
        
        # A mapping from neighborhood patterns to their index in the rule's binary string.
        # e.g., '111' corresponds to index 0, '110' to index 1, etc., following the Wolfram convention.
        neighborhoods = {
            (1, 1, 1): 0, (1, 1, 0): 1, (1, 0, 1): 2, (1, 0, 0): 3,
            (0, 1, 1): 4, (0, 1, 0): 5, (0, 0, 1): 6, (0, 0, 0): 7
        }
        
        state = [int(c) for c in state_str]
        length = len(state)
        next_state = [0] * length

        for i in range(length):
            # Determine the 3-cell neighborhood, wrapping around the edges.
            left = state[i - 1]
            center = state[i]
            right = state[(i + 1) % length]
            
            pattern = (left, center, right)
            
            # Find the output for this pattern from the rule.
            rule_index = neighborhoods[pattern]
            next_state[i] = int(rule_bin[rule_index])
            
        return "".join(map(str, next_state))

    valid_intermediate_steps = set()
    
    # Iterate through all 256 possible elementary cellular automaton rules.
    for rule in range(256):
        # First evolution: gen0 -> gen1
        gen1_candidate = evolve(gen0_str, rule)
        
        # Second evolution: gen1 -> gen2
        gen2_candidate = evolve(gen1_candidate, rule)
        
        # Check if the result matches the given final generation
        if gen2_candidate == gen2_str:
            valid_intermediate_steps.add(gen1_candidate)

    # The problem implies a unique solution, so we expect to find exactly one.
    if len(valid_intermediate_steps) == 1:
        solution = valid_intermediate_steps.pop()
        # Output the single valid solution for the intermediate step.
        print(solution)
    elif len(valid_intermediate_steps) == 0:
        print("No solution found.")
    else:
        # This case is not expected based on the problem description.
        print("Multiple possible intermediate steps found:")
        for s in sorted(list(valid_intermediate_steps)):
            print(s)

solve_automaton()