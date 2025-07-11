import collections

def has_glider_for_rule(rule):
    """
    Checks if a given ECA rule supports any gliders by simulating the evolution of various small initial patterns.
    """
    max_initial_width = 8  # Maximum width of initial patterns to test
    max_steps = 400        # Maximum simulation steps for each pattern
    max_pattern_width = 150 # Stop if pattern grows too large (likely chaotic)

    # --- Helper function to perform one step of ECA evolution ---
    def evolve(config, rule_code):
        padded_config = [0, 0] + config + [0, 0]
        next_config = []
        for i in range(len(padded_config) - 2):
            neighborhood = tuple(padded_config[i:i+3])
            index = 4 * neighborhood[0] + 2 * neighborhood[1] + 1 * neighborhood[2]
            output = (rule_code >> index) & 1
            next_config.append(output)
        return next_config

    # --- Helper function to trim zero-padding and find pattern's position ---
    def trim_config(config):
        try:
            first_one = config.index(1)
            last_one = len(config) - 1 - config[::-1].index(1)
            pattern = config[first_one:last_one + 1]
            position = first_one
            return pattern, position
        except ValueError:
            return [], 0 # All zeros

    # --- Main loop to test different initial patterns ---
    for width in range(1, max_initial_width + 1):
        # Iterate through all non-trivial binary patterns of the given width
        for i in range(1, 1 << width):
            initial_pattern = [int(b) for b in bin(i)[2:].zfill(width)]
            
            current_pattern = initial_pattern
            current_pos = 0
            
            # History maps a pattern shape to the (step, position) it was first seen at.
            # This helps detect cycles and gliders.
            history = {}

            for step in range(max_steps):
                pattern_tuple = tuple(current_pattern)
                
                # Check if this pattern has been seen before
                if pattern_tuple in history:
                    prev_step, prev_pos = history[pattern_tuple]
                    period = step - prev_step
                    displacement = current_pos - prev_pos
                    if displacement != 0 and period > 0:
                        # Found a configuration that repeats at a different location. It's a glider.
                        return True
                    else:
                        # It's a stable or blinking pattern, not a glider. Stop this simulation.
                        break
                
                history[pattern_tuple] = (step, current_pos)

                # Evolve to the next step
                next_config_untrimmed = evolve(current_pattern, rule)
                
                # Process the new configuration
                next_pattern, trim_offset = trim_config(next_config_untrimmed)
                
                # Check for termination conditions
                if not next_pattern:  # Pattern died out
                    break
                if len(next_pattern) > max_pattern_width:  # Pattern growing too large
                    break
                
                # Calculate the position of the new pattern. The evolution itself causes a shift of -1.
                next_pos = current_pos - 1 + trim_offset
                
                # Update state for the next iteration
                current_pattern = next_pattern
                current_pos = next_pos
                
    return False # No glider found for this rule with the tested patterns

def solve():
    """
    Finds and counts all compact ECAs that have at least one glider.
    """
    glider_rules = []
    # An ECA is compact if the all-zeroes neighborhood maps to 0.
    # This means its rule number must be even.
    for rule in range(0, 256, 2):
        if has_glider_for_rule(rule):
            glider_rules.append(rule)
    
    print("The following compact Elementary Cellular Automata have gliders:")
    # "Remember in the final code you still need to output each number in the final equation!"
    # I will interpret this as printing the list of rule numbers found.
    print(glider_rules)
    print("\nTotal count:")
    print(len(glider_rules))

if __name__ == '__main__':
    solve()