import collections

def has_glider(rule_num):
    """
    Searches for a glider in a given ECA rule.

    This function tests a set of initial small patterns, evolves them,
    and checks if any pattern repeats itself at a shifted position.
    """
    # Search parameters
    MAX_INITIAL_WIDTH = 7  # Test initial patterns up to this width
    MAX_STEPS = 256        # Simulate for this many steps
    MAX_PATTERN_WIDTH = 100 # Stop if pattern grows too large

    # --- Step 1: Set up the rule mapping ---
    rule_map = {}
    binary_rule = format(rule_num, '08b')
    for i in range(8):
        # Neighborhoods are ordered from '111' (7) down to '000' (0)
        neighborhood = tuple(int(bit) for bit in format(7 - i, '03b'))
        rule_map[neighborhood] = int(binary_rule[i])

    # --- Step 2: Iterate through initial seed patterns ---
    # We test 2^MAX_INITIAL_WIDTH - 1 different starting patterns.
    for i in range(1, 2**MAX_INITIAL_WIDTH):
        pattern_str = bin(i)[2:]
        config = tuple(int(b) for b in pattern_str)
        
        # History tracks {pattern: [(time, position), ...]}
        history = collections.defaultdict(list)
        position = 0

        # --- Step 3: Simulate the evolution of the seed ---
        for t in range(MAX_STEPS):
            # Check for termination conditions
            if not config or len(config) > MAX_PATTERN_WIDTH:
                break
            
            # Record current state in history
            history[config].append((t, position))

            # Evolve the configuration
            # Pad the config to calculate the next state of boundary cells
            padded_config = (0, 0) + config + (0, 0)
            next_config_list = []
            for j in range(len(padded_config) - 2):
                neighborhood = padded_config[j:j+3]
                next_config_list.append(rule_map[neighborhood])

            # Trim leading/trailing zeros to get the new pattern and its shift
            try:
                first_one = next_config_list.index(1)
                last_one = len(next_config_list) - 1 - next_config_list[::-1].index(1)
                
                # The shift is relative to the previous pattern's left edge
                shift = first_one - 2
                position += shift
                config = tuple(next_config_list[first_one : last_one + 1])
            except ValueError:
                # The pattern has vanished (all zeros)
                config = tuple()

            # --- Step 4: Check for glider condition ---
            # If the new pattern has been seen before
            if config in history:
                for prev_t, prev_pos in history[config]:
                    # If position is different, it's a glider
                    if position != prev_pos:
                        return True
                        
    return False

def find_compact_ecas_with_gliders():
    """
    Main function to find and count compact ECAs with gliders.
    """
    print("Searching for compact ECAs with gliders...")
    print("This may take several minutes.")
    
    glider_rules = []
    # Iterate through all 256 ECA rules
    for rule_num in range(256):
        # A compact ECA must be an even-numbered rule (output for '000' is 0)
        if rule_num % 2 == 0:
            if has_glider(rule_num):
                glider_rules.append(rule_num)
    
    print("\nFound the following compact ECAs with gliders:")
    print(glider_rules)
    print(f"\nTotal number of compact ECAs with a glider: {len(glider_rules)}")
    return len(glider_rules)

if __name__ == '__main__':
    count = find_compact_ecas_with_gliders()
    # The final answer format
    print(f"\n<<<{count}>>>")
