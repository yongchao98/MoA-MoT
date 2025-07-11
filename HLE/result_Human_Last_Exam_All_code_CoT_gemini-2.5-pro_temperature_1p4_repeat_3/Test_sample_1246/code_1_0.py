def get_rule_map(rule_num):
    """Generates the rule mapping for a given ECA rule number."""
    b_string = format(rule_num, '08b')
    rule_map = {}
    for i in range(8):
        neighborhood_val = 7 - i
        neighborhood = tuple(int(x) for x in format(neighborhood_val, '03b'))
        output = int(b_string[i])
        rule_map[neighborhood] = output
    return rule_map

def evolve(config, rule_map):
    """Evolves a compact configuration one step forward."""
    # Pad with two zeros on each side for boundary conditions.
    padded = [0, 0] + config + [0, 0]
    next_config = []
    for i in range(len(padded) - 2):
        neighborhood = tuple(padded[i:i+3])
        next_config.append(rule_map.get(neighborhood, 0))
    return next_config

def trim_config(config):
    """Removes leading/trailing zeros from a configuration and returns the number of leading zeros removed."""
    if not any(config):
        return [], 0
    
    start_index = -1
    end_index = -1
    for i, val in enumerate(config):
        if val == 1:
            if start_index == -1:
                start_index = i
            end_index = i
            
    if start_index == -1:
        return [], 0
        
    return config[start_index:end_index+1], start_index

def has_glider(rule_num, max_width, max_time):
    """
    Searches for a glider in a given compact ECA rule. This is a heuristic search.
    """
    rule_map = get_rule_map(rule_num)

    # Iterate through a set of simple, non-trivial initial configurations.
    for width in range(1, max_width + 1):
        # Generate all 2^(w-2) configurations of 'width' that start and end with '1'.
        num_configs_for_width = 2**max(0, width - 2)
        
        for i in range(num_configs_for_width):
            if width == 1:
                initial_core = [1]
            else:
                middle_pattern = [int(x) for x in format(i, f'0{width-2}b')]
                initial_core = [1] + middle_pattern + [1]
            
            # --- Simulation loop for this initial configuration ---
            current_config = initial_core
            current_pos = 0
            
            # History tracks seen shapes and their (time, position)
            history = {tuple(initial_core): [(0, 0)]}

            for t in range(1, max_time + 1):
                untrimmed_next = evolve(current_config, rule_map)
                
                if not any(untrimmed_next): # Configuration died out
                    break
                
                next_core, leading_zeros = trim_config(untrimmed_next)
                
                # Evolution shifts frame left by 1; trimming shifts it right.
                next_pos = current_pos - 1 + leading_zeros
                
                core_tuple = tuple(next_core)
                if core_tuple in history:
                    for _old_t, old_pos in history[core_tuple]:
                        # If shape reappears at a new position, it's a glider.
                        if next_pos != old_pos:
                            return True
                
                if core_tuple not in history:
                    history[core_tuple] = []
                history[core_tuple].append((t, next_pos))

                current_config = next_core
                current_pos = next_pos
                
                # Heuristic break for overly complex or large patterns to manage runtime
                if len(history) > max_time * 2 or len(current_config) > max_width * 5:
                     break
    return False

def solve():
    """
    Finds the number of compact Elementary Cellular Automata that have a glider.
    """
    # The search for gliders is computationally intensive. These parameters control
    # the search depth. They are a compromise for finding many known gliders in
    # a reasonable time.
    MAX_INITIAL_WIDTH = 7
    MAX_SIMULATION_TIME = 128

    glider_rules = []
    for rule in range(256):
        # A "compact" ECA requires that the rule for (0,0,0) is 0.
        # This corresponds to the rule number being even.
        if rule % 2 == 0:
            if has_glider(rule, MAX_INITIAL_WIDTH, MAX_SIMULATION_TIME):
                glider_rules.append(rule)

    print(f"Based on the search, {len(glider_rules)} compact ECAs have a glider.")
    print("The final count is determined by the length of the list of rules found.")
    for rule in glider_rules:
        print(f"Rule {rule} was found to have a glider.")
    
    print(f"\nFinal Answer: {len(glider_rules)}")


if __name__ == '__main__':
    solve()
<<<16>>>