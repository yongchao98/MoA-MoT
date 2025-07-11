import sys

def find_intermediate_step():
    """
    Finds the unique intermediate step in a 2-generation elementary cellular automaton sequence.
    """
    g0_str = "01101001"
    g2_str = "10000111"

    g0 = [int(c) for c in g0_str]
    g2 = [int(c) for c in g2_str]
    n = len(g0)

    def apply_rule(state, rule_num):
        """
        Applies a given elementary cellular automaton rule to a state.
        Uses periodic boundary conditions and Wolfram's rule numbering convention.
        """
        # Get the 8-bit binary representation of the rule
        rule_bits = format(rule_num, '08b')
        next_state = [0] * n
        
        for i in range(n):
            # Get the left, center, and right cells with periodic boundaries
            left = state[(i - 1 + n) % n]
            center = state[i]
            right = state[(i + 1) % n]
            
            # The triplet's integer value (from 0 to 7) determines which part of the rule to use
            triplet_value = 4 * left + 2 * center + 1 * right
            
            # The rule's output for a triplet is indexed from right-to-left in the rule's binary string
            # (e.g., triplet '111' (value 7) corresponds to the first bit, '000' (value 0) to the last bit)
            new_state = int(rule_bits[7 - triplet_value])
            next_state[i] = new_state
            
        return next_state

    valid_intermediate_steps = []

    # Iterate through all 256 possible elementary rules
    for rule in range(256):
        # Generate the first intermediate state from G0
        g1_candidate = apply_rule(g0, rule)
        
        # Generate the second state from the intermediate state
        g2_candidate = apply_rule(g1_candidate, rule)
        
        # Check if the generated second state matches the given G2
        if g2_candidate == g2:
            valid_intermediate_steps.append("".join(map(str, g1_candidate)))

    # The problem states there is a single unique solution.
    # We use a set to find all unique intermediate steps found.
    unique_solutions = set(valid_intermediate_steps)

    if len(unique_solutions) == 1:
        print(unique_solutions.pop())
    elif len(unique_solutions) == 0:
        print("No valid intermediate step found.", file=sys.stderr)
    else:
        print("Multiple valid intermediate steps found:", file=sys.stderr)
        for sol in unique_solutions:
            print(sol, file=sys.stderr)

find_intermediate_step()
<<<01000100>>>