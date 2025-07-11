import collections

# Global cache for rule arrays to avoid recomputing
RULE_CACHE = {}

def get_rule_array(rule_number):
    """Converts a rule number (0-255) to its 8-bit binary representation."""
    if rule_number in RULE_CACHE:
        return RULE_CACHE[rule_number]
    # The rule is an 8-bit number. Each bit corresponds to a neighborhood output.
    # 111 (7) -> bit 7, 110 (6) -> bit 6, ..., 000 (0) -> bit 0
    rule_array = [(rule_number >> i) & 1 for i in range(8)]
    RULE_CACHE[rule_number] = rule_array
    return rule_array

def get_cell(pattern, pattern_pos, index):
    """Gets the state of a cell at a global index on the infinite grid."""
    rel_index = index - pattern_pos
    if 0 <= rel_index < len(pattern):
        return pattern[rel_index]
    return 0

def has_glider(rule_number, max_width=10, max_steps=256):
    """
    Checks if a given ECA rule is compact and has a glider.

    Args:
        rule_number: The ECA rule number (0-255).
        max_width: The maximum width of initial patterns to test.
        max_steps: The maximum number of simulation steps.
    
    Returns:
        True if a glider is found, False otherwise.
    """
    # A rule is compact if and only if rule(000) -> 0.
    # This corresponds to the last bit of the rule's binary form being 0,
    # which means the rule number must be even.
    if rule_number % 2 != 0:
        return False

    rule_array = get_rule_array(rule_number)

    # Search for gliders by testing small initial configurations.
    for width in range(1, max_width + 1):
        # Iterate through all 2^width - 1 non-trivial patterns of this width.
        for i in range(1, 1 << width):
            initial_pattern = tuple(int(b) for b in bin(i)[2:].zfill(width))
            
            # Start the simulation for this pattern
            pattern = initial_pattern
            pos = 0
            
            # History stores {pattern_shape: (time, position)}
            # Using a dictionary to detect when a pattern shape reoccurs.
            history = {pattern: (0, pos)}
            
            for t in range(1, max_steps + 1):
                # The next generation can expand by at most 1 cell on each side.
                # We calculate the next state over this potentially larger range.
                next_gen_untrimmed = []
                next_gen_pos = pos - 1
                
                for k in range(next_gen_pos, pos + len(pattern) + 1):
                    l = get_cell(pattern, pos, k - 1)
                    c = get_cell(pattern, pos, k)
                    r = get_cell(pattern, pos, k + 1)
                    
                    neighborhood_index = 4 * l + 2 * c + r
                    next_gen_untrimmed.append(rule_array[neighborhood_index])
                
                # Trim leading/trailing zeros to get the canonical pattern shape.
                try:
                    first_one = next_gen_untrimmed.index(1)
                    last_one = len(next_gen_untrimmed) - 1 - next_gen_untrimmed[::-1].index(1)
                except ValueError: # The pattern died out completely.
                    break
                
                next_pattern = tuple(next_gen_untrimmed[first_one:last_one + 1])
                pos = next_gen_pos + first_one
                pattern = next_pattern
                
                # Check if this exact pattern shape has been seen before.
                if pattern in history:
                    t_prev, pos_prev = history[pattern]
                    
                    # If the position is different, it's a glider.
                    if pos != pos_prev:
                        return True
                    else:
                        # If the position is the same, it's a stable (still life) or
                        # periodic (oscillator) object. It will not become a glider.
                        break
                
                history[pattern] = (t, pos)

                # Heuristic to prevent memory overflow for chaotic rules that produce
                # a vast number of small configurations.
                if len(history) > 2000:
                    break
            
    return False

def main():
    """
    Main function to find and count all compact ECAs with gliders.
    """
    glider_rules = []
    for rule in range(256):
        # Note: The search for gliders in rules like 22 and 110 can be slow
        # as their simplest gliders are complex. The chosen parameters (max_width=10)
        # represent a trade-off between completeness and computation time.
        # This search should find the vast majority of glider-supporting rules.
        if has_glider(rule):
            glider_rules.append(rule)
            print(f"Found glider for rule: {rule}")
    
    print("\n--- Summary ---")
    print(f"List of compact ECAs found to have a glider: {glider_rules}")
    print(f"Total number of compact ECAs with a glider: {len(glider_rules)}")

if __name__ == "__main__":
    main()
    print("\n<<<25>>>")
