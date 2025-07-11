import sys

def find_compact_ecas_with_gliders():
    """
    Finds the number of compact Elementary Cellular Automata that have a glider.

    This function iterates through all 256 ECA rules, identifies the compact ones,
    and then searches for gliders by simulating small initial configurations.
    """
    compact_ecas_with_gliders = set()

    # --- Search Parameters ---
    # MAX_INITIAL_WIDTH: The maximum width of the initial pattern to test.
    # MAX_STEPS: The maximum number of simulation steps to check for a glider.
    # These parameters are chosen to be large enough to find common gliders
    # while keeping the computation time reasonable.
    MAX_INITIAL_WIDTH = 8
    MAX_STEPS = 128

    # Pre-calculate rule maps for all 256 rules for efficiency.
    rule_maps = []
    for rule_num in range(256):
        rule_map = {}
        binary_rule = format(rule_num, '08b')
        for i in range(8):
            neighborhood_int = 7 - i
            neighborhood_tuple = tuple(int(b) for b in format(neighborhood_int, '03b'))
            rule_map[neighborhood_tuple] = int(binary_rule[i])
        rule_maps.append(rule_map)

    # Iterate through all 256 ECA rules.
    for rule_num in range(256):
        # A rule is compact if and only if it maps the '000' neighborhood to '0'.
        # This corresponds to even rule numbers.
        if rule_num % 2 != 0:
            continue
        
        # Rule 0 annihilates all patterns, so it cannot have a non-trivial glider.
        if rule_num == 0:
            continue

        rule_map = rule_maps[rule_num]
        found_glider_for_this_rule = False

        # Iterate through initial configurations of different widths.
        for width in range(1, MAX_INITIAL_WIDTH + 1):
            # Iterate through all non-trivial patterns of the given width.
            for i in range(1, 2**width):
                initial_pattern_list = [int(b) for b in bin(i)[2:].zfill(width)]
                initial_pattern_tuple = tuple(initial_pattern_list)
                
                # Set up the simulation grid with padding.
                padding = MAX_STEPS
                grid_size = width + 2 * padding
                grid = [0] * grid_size
                initial_start_pos = padding
                grid[initial_start_pos : initial_start_pos + width] = initial_pattern_list
                
                current_grid = grid
                
                # History to detect non-glider cycles (still lifes/oscillators).
                history = {(initial_pattern_tuple, initial_start_pos)}

                # Simulate the evolution for MAX_STEPS.
                for step in range(1, MAX_STEPS + 1):
                    next_grid = [0] * grid_size
                    
                    # Find pattern bounds to optimize the simulation step.
                    first_one = -1
                    last_one = -1
                    try:
                        first_one = current_grid.index(1)
                        # Find last '1' by searching from the end.
                        last_one = grid_size - 1 - current_grid[::-1].index(1)
                    except ValueError: # No '1's left in the grid
                        break

                    # Evolve the relevant part of the grid.
                    for j in range(first_one - 1, last_one + 2):
                        if 1 <= j < grid_size - 1:
                            neighborhood = tuple(current_grid[j-1 : j+2])
                            next_grid[j] = rule_map[neighborhood]
                    
                    current_grid = next_grid
                    
                    # Extract the new pattern to check for gliders.
                    new_first_one = -1
                    try:
                        new_first_one = current_grid.index(1)
                        new_last_one = grid_size - 1 - current_grid[::-1].index(1)
                    except ValueError: # Pattern disappeared
                        break

                    current_pattern_tuple = tuple(current_grid[new_first_one : new_last_one + 1])
                    
                    # Check if it's a glider (same pattern, different location).
                    if current_pattern_tuple == initial_pattern_tuple:
                        displacement = new_first_one - initial_start_pos
                        if displacement != 0:
                            compact_ecas_with_gliders.add(rule_num)
                            found_glider_for_this_rule = True
                            break
                    
                    # Prune search by detecting cycles.
                    current_state = (current_pattern_tuple, new_first_one)
                    if current_state in history:
                        break 
                    history.add(current_state)

                if found_glider_for_this_rule:
                    break
            if found_glider_for_this_rule:
                break
    
    # Print the results
    sorted_rules = sorted(list(compact_ecas_with_gliders))
    print(f"Found {len(sorted_rules)} compact ECAs with gliders.")
    print("The rule numbers are:")
    # The problem asks to output each number
    for rule in sorted_rules:
        print(rule, end=' ')
    print()

find_compact_ecas_with_gliders()
<<<45>>>