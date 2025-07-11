import sys

def solve_eca_gliders():
    """
    This script finds the number of compact Elementary Cellular Automata (ECA)
    that have at least one glider.

    - A 'compact' ECA is one that keeps configurations with a finite number of 1s finite.
      This corresponds to all even-numbered rules.
    - A 'glider' is a finite pattern that repeats itself at a new location after some time.

    The script iterates through all even rules and simulates a variety of small starting
    patterns to detect glider-like behavior.
    """

    # Search parameters defining the scope of our glider search.
    # We will test all initial patterns of width up to MAX_INIT_WIDTH.
    MAX_INIT_WIDTH = 12
    # We will simulate each pattern for up to MAX_STEPS.
    MAX_STEPS = 512
    # To avoid infinite loops on chaotic rules, we stop if a pattern's width exceeds this limit.
    MAX_PATTERN_WIDTH = 100

    memo = {} # Memoization to cache results for each rule.

    def has_glider(rule_num):
        """
        Checks if a given ECA rule has a glider by simulating various initial patterns.
        """
        if rule_num in memo:
            return memo[rule_num]

        # A rule must be even-numbered to be compact.
        if rule_num % 2 != 0:
            memo[rule_num] = False
            return False

        # Precompute the rule's transition table for fast lookups.
        rule_map = {}
        for i in range(8):
            neighborhood = tuple(int(b) for b in format(i, '03b'))
            rule_map[neighborhood] = (rule_num >> i) & 1

        # Iterate through a set of initial configurations (patterns).
        # Integers from 1 up to 2**MAX_INIT_WIDTH are used, where their binary
        # representation defines the starting pattern of 1s and 0s.
        for i in range(1, 1 << MAX_INIT_WIDTH):
            initial_pattern = tuple(int(b) for b in bin(i)[2:])
            
            current_pattern = initial_pattern
            current_pos = 0  # Position of the pattern's first cell
            
            # History tracks patterns seen during this simulation.
            # Format: {pattern_tuple: (time_step, position)}
            history = {current_pattern: (0, current_pos)}

            for t in range(1, MAX_STEPS):
                # Evolve the pattern for one time step.
                # Pad the pattern to simulate an infinite grid of 0s.
                padded = (0, 0) + current_pattern + (0, 0)
                next_gen_list = [rule_map[tuple(padded[j:j+3])] for j in range(len(padded) - 2)]

                # Process the new generation to find the pattern and its position.
                try:
                    first_one = next_gen_list.index(1)
                    last_one = len(next_gen_list) - 1 - next_gen_list[::-1].index(1)
                except ValueError:
                    # The pattern vanished (all 0s).
                    break
                
                next_pattern = tuple(next_gen_list[first_one:last_one+1])
                
                # The position of the new pattern is relative to the old one.
                # The evolution shifts the coordinate system by -1.
                next_pos = current_pos - 1 + first_one

                # Heuristic to stop simulating chaotically growing patterns.
                if len(next_pattern) > MAX_PATTERN_WIDTH:
                    break

                # Check history to see if this exact pattern has appeared before.
                if next_pattern in history:
                    t_prev, pos_prev = history[next_pattern]
                    
                    displacement = next_pos - pos_prev
                    if displacement != 0:
                        # Glider found: the pattern repeated at a different location.
                        memo[rule_num] = True
                        return True
                
                # Store the new pattern in history and continue.
                history[next_pattern] = (t, next_pos)
                current_pattern = next_pattern
                current_pos = next_pos
        
        # No glider found within the search limits for this rule.
        memo[rule_num] = False
        return False

    glider_rules = []
    for rule_num in range(256):
        if has_glider(rule_num):
            glider_rules.append(rule_num)

    print("The following compact ECA rules were found to have a glider:")
    print(glider_rules)
    print("\nThe problem asks for the total count of such rules. The numbers of the rules found are shown above.")
    print(f"\nThe final count is {len(glider_rules)}.")

    return len(glider_rules)

# Run the solver and capture the final numerical answer.
final_count = solve_eca_gliders()
print(f"\n<<<>>>\n{final_count}")