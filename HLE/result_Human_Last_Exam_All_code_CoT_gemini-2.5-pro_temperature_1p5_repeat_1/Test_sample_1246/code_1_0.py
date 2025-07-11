def find_compact_ecas_with_gliders():
    """
    This script analyzes all 256 Elementary Cellular Automata (ECA) to identify
    which ones are "compact" and have a "glider" as per the specified definitions.

    - A "compact" ECA must map the '000' neighborhood to '0', meaning its
      rule number must be even.
    - A "glider" is a finite pattern that repeats itself at a new location after
      a number of steps.

    The script iterates through all even-numbered rules, and for each rule, it
    simulates a variety of small initial configurations to detect glider-like behavior.
    """

    # --- Simulation Parameters ---
    # Test initial patterns represented by integers up to this value.
    # e.g., 255 tests all non-trivial patterns of width 1 to 8.
    MAX_INIT_CONFIG_INT = 255
    # Maximum number of simulation steps for any single starting pattern.
    MAX_STEPS = 512
    # If a pattern's width exceeds this, we stop, assuming it's not a stable glider.
    MAX_WIDTH = 64
    # Padding applied to the configuration edges to correctly handle neighborhoods.
    PADDING = 2

    compact_glider_rules = []

    # 1. Iterate through all 256 ECA rules.
    for rule_num in range(256):
        # A rule is compact only if it's even (000 -> 0).
        if rule_num % 2 != 0:
            continue

        # Generate the rule's transition map (neighborhood -> output).
        rule_bits = format(rule_num, '08b')
        rule_map = {
            tuple(int(b) for b in format(7 - i, '03b')): int(rule_bits[i])
            for i in range(8)
        }

        found_glider_for_rule = False

        # 2. Test a range of small initial configurations for this rule.
        for i in range(1, MAX_INIT_CONFIG_INT + 1):
            if found_glider_for_rule:
                break

            # Initial state: a compact configuration and its position.
            config = [int(b) for b in bin(i)[2:]]
            position = 0
            history = {}  # Stores {pattern_tuple: (time, position)}

            # 3. Simulate the evolution.
            for t in range(MAX_STEPS):
                if not config or len(config) > MAX_WIDTH:
                    break  # Pattern vanished or grew too large.

                config_tuple = tuple(config)

                # 4. Check for a repeating pattern.
                if config_tuple in history:
                    prev_t, prev_pos = history[config_tuple]
                    shift = position - prev_pos

                    # A glider must have moved (shift is non-zero).
                    if shift != 0:
                        compact_glider_rules.append(rule_num)
                        found_glider_for_rule = True
                        break
                
                history[config_tuple] = (t, position)

                # Evolve the configuration by one step.
                padded_config = ([0] * PADDING) + config + ([0] * PADDING)
                next_config_padded = [
                    rule_map[tuple(padded_config[j:j+3])]
                    for j in range(len(padded_config) - 2)
                ]
                
                # Trim leading/trailing zeros to keep the config compact.
                try:
                    first_one = next_config_padded.index(1)
                    last_one = len(next_config_padded) - 1 - next_config_padded[::-1].index(1)
                    
                    # Update config and position for the next step.
                    position += first_one - PADDING
                    config = next_config_padded[first_one : last_one + 1]
                except ValueError:
                    config = [] # No '1's found, pattern vanished.
    
    return sorted(compact_glider_rules)


if __name__ == '__main__':
    glider_rules = find_compact_ecas_with_gliders()
    count = len(glider_rules)
    
    print(f"Found {count} compact ECAs that have a glider.")
    print("The rule numbers are:")
    # Printing the list of numbers in the requested "equation" format.
    print(" + ".join(map(str, glider_rules)))
    print(f"\nThe total count is {count}.")
