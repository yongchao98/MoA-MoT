import sys

# It is recommended to run this script using PyPy for a significant speed-up.

def solve():
    """
    Finds and counts the number of compact Elementary Cellular Automata (ECAs) that have a glider.

    This script performs a computational search based on the problem's definitions:
    - A 'compact ECA' is an even-numbered rule, ensuring a 0-background remains a 0-background.
    - A 'glider' is a pattern that repeats itself at a new location after some time.

    The search is computationally intensive. The parameters below are chosen based on published
    research to be sufficient for finding known gliders.
    - max_width: The maximum width of initial patterns (seeds) to test.
    - max_steps: The maximum number of simulation steps for each seed.
    - max_pattern_len: A limit to stop simulating patterns that grow too large (likely chaotic).

    An optimization is used: if rule R has a glider, its reflection R' does too,
    halving the search space.
    """

    MAX_WIDTH = 12
    MAX_STEPS = 1024
    MAX_PATTERN_LEN = 100

    # Cache for rule binary representations
    RULE_BITS_CACHE = {i: format(i, '08b') for i in range(256)}

    # Cache for reflection mapping to avoid re-computation
    REFLECTION_CACHE = {}

    def get_reflected_rule(rule_num):
        """Calculates the rule number for the mirror image of a given rule."""
        if rule_num in REFLECTION_CACHE:
            return REFLECTION_CACHE[rule_num]

        b = RULE_BITS_CACHE[rule_num]
        # Standard reflection mapping for ECA rules
        # new_bit_k = old_bit_rev(k)
        # e.g., new_bit_6 (for 110) = old_bit_3 (for 011)
        # In our b7...b0 string format, bit k corresponds to index (7-k)
        # However, the source bits are usually referred to by their value.
        # k=0->0, k=1->4, k=2->2, k=3->6, k=4->1, k=5->5, k=6->3, k=7->7
        bits_map = [b[0], b[4], b[2], b[6], b[1], b[5], b[3], b[7]]
        reflected_num = int("".join(bits_map), 2)
        REFLECTION_CACHE[rule_num] = reflected_num
        return reflected_num


    def find_glider(rule_num):
        """
        Searches for a glider in the given ECA rule.

        Returns True if a glider is found, False otherwise.
        """
        rule_bits = RULE_BITS_CACHE[rule_num]
        rule_map = {
            (1, 1, 1): int(rule_bits[0]), (1, 1, 0): int(rule_bits[1]),
            (1, 0, 1): int(rule_bits[2]), (1, 0, 0): int(rule_bits[3]),
            (0, 1, 1): int(rule_bits[4]), (0, 1, 0): int(rule_bits[5]),
            (0, 0, 1): int(rule_bits[6]), (0, 0, 0): int(rule_bits[7]),
        }

        # Iterate through all possible non-trivial initial seeds up to MAX_WIDTH
        for i in range(1, 1 << MAX_WIDTH):
            initial_pattern = tuple(int(b) for b in bin(i)[2:])
            
            # Start simulation for this seed
            current_config = initial_pattern
            current_offset = 0
            
            # History stores pattern -> (time, offset)
            history = {current_config: (0, current_offset)}

            for t in range(1, MAX_STEPS + 1):
                config_list = list(current_config)
                # Pad the grid to correctly calculate the next state of boundary cells
                grid = [0, 0] + config_list + [0, 0]
                
                next_raw = [rule_map[tuple(grid[j-1:j+2])] for j in range(1, len(grid)-1)]
                
                # Stop if pattern vanishes
                if 1 not in next_raw:
                    break

                # Trim leading/trailing zeros for a canonical representation
                first_one = next_raw.index(1)
                last_one = len(next_raw) - 1 - next_raw[::-1].index(1)
                
                next_config = tuple(next_raw[first_one : last_one + 1])
                
                # The grid for the next step is shifted left by 1 relative to the old one.
                # The first '1' in the new raw config gives its position within this new grid.
                next_offset = current_offset - 1 + first_one

                # Check if this exact pattern has been seen before
                if next_config in history:
                    prev_t, prev_offset = history[next_config]
                    shift = next_offset - prev_offset
                    
                    # If it's at a different location, it's a glider
                    if shift != 0:
                        return True
                    # If it's at the same location, it's a stable or oscillating pattern.
                    # No need to simulate this seed further.
                    else:
                        break

                # Prune search for patterns that grow too large
                if len(next_config) > MAX_PATTERN_LEN:
                    break
                
                # Record the new pattern and continue
                history[next_config] = (t, next_offset)
                current_config = next_config
                current_offset = next_offset
        
        return False

    glider_rules = set()
    checked_rules = set()

    # Iterate through all compact ECAs (even-numbered rules)
    for rule_num in range(0, 256, 2):
        if rule_num in checked_rules:
            continue

        # Check for progress
        # sys.stdout.write(f"\rChecking rule {rule_num}...")
        # sys.stdout.flush()

        if find_glider(rule_num):
            glider_rules.add(rule_num)
            # Its reflection will also have gliders
            reflected = get_reflected_rule(rule_num)
            glider_rules.add(reflected)
        
        # Mark both the rule and its reflection as checked
        checked_rules.add(rule_num)
        checked_rules.add(get_reflected_rule(rule_num))

    # sys.stdout.write("\nDone.\n")
    print(len(glider_rules))


solve()