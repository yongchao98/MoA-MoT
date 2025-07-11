import sys

def find_compact_ecas_with_gliders():
    """
    This script finds the number of compact Elementary Cellular Automata (ECAs) that have a glider.
    - An ECA is 'compact' if its rule number is even (meaning the '000' neighborhood maps to 0).
    - A 'glider' is a finite pattern that repeats its shape at a new location after some time.
    - The script iterates through all 128 compact ECAs and searches for a glider by simulating
      small initial patterns.
    """

    # --- Search Parameters ---
    # These parameters are a trade-off between search completeness and computation time.
    # Higher values will find more complex gliders but will take much longer to run.
    MAX_INIT_WIDTH = 7  # Maximum width of initial patterns to test.
    MAX_STEPS = 70      # Maximum simulation steps to check for a repeating pattern.
    # -------------------------
    
    # A grid wide enough to prevent boundary effects during simulation.
    GRID_WIDTH = MAX_INIT_WIDTH + 2 * MAX_STEPS + 4

    def get_rule_output(rule, neighborhood_index):
        """Gets the output of a rule for a given neighborhood index."""
        return (rule >> neighborhood_index) & 1

    def evolve_step(config, rule_outputs):
        """Performs one step of the ECA evolution on a configuration."""
        new_config = [0] * len(config)
        for i in range(1, len(config) - 1):
            neighborhood_index = 4 * config[i - 1] + 2 * config[i] + 1 * config[i + 1]
            new_config[i] = rule_outputs[neighborhood_index]
        return new_config

    def trim_config(config):
        """Removes leading/trailing zeros and returns the trimmed pattern and its start position."""
        try:
            first_one = config.index(1)
            last_one = len(config) - 1 - config[::-1].index(1)
            return config[first_one : last_one + 1], first_one
        except ValueError:  # All zeros
            return [], 0

    def has_glider(rule):
        """Searches for a glider for a given rule within the specified search limits."""
        # Pre-calculate rule outputs for speed
        rule_outputs = [get_rule_output(rule, i) for i in range(8)]
        
        # Iterate through initial patterns (represented by integers 1 to 2^width - 1)
        for i in range(1, 2**MAX_INIT_WIDTH):
            init_pattern = [int(x) for x in bin(i)[2:]]
            
            # Setup the grid with the pattern in the center
            grid = [0] * GRID_WIDTH
            start_pos = (GRID_WIDTH - len(init_pattern)) // 2
            grid[start_pos : start_pos + len(init_pattern)] = init_pattern

            # Simulate for MAX_STEPS
            current_grid = grid
            for _ in range(MAX_STEPS):
                current_grid = evolve_step(current_grid, rule_outputs)
                trimmed_pattern, new_pos = trim_config(current_grid)

                if not trimmed_pattern:
                    break  # Pattern died out

                # Check if the shape matches the initial pattern
                if trimmed_pattern == init_pattern:
                    shift = new_pos - start_pos
                    if shift != 0:
                        return True  # Glider found
        return False # No glider found within search limits

    compact_glider_rules = []
    # Iterate through all even rules (compact ECAs) from 0 to 254.
    for rule in range(0, 256, 2):
        sys.stdout.write(f"\rChecking Rule {rule}...")
        sys.stdout.flush()
        if has_glider(rule):
            compact_glider_rules.append(rule)

    print("\r" + " " * 30) # Clear the line
    print("Search complete.")
    print("\nThe compact ECA rules found to have a glider are:")
    # The prompt requires printing each number
    print(", ".join(map(str, compact_glider_rules)))
    
    count = len(compact_glider_rules)
    print(f"\nTotal number of compact ECAs with a glider found: {count}")
    return count

if __name__ == '__main__':
    final_count = find_compact_ecas_with_gliders()
    print(f"\n<<< {final_count} >>>")