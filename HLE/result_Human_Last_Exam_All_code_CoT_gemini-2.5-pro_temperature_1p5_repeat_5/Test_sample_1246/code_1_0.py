import sys

def solve():
    """
    Finds the number of compact Elementary Cellular Automata that have a glider.
    """
    # Configuration for the search. These bounds are a trade-off between
    # completeness and computational time.
    MAX_INIT_WIDTH = 8  # Maximum width of initial patterns to test.
    MAX_STEPS = 128     # Maximum number of steps to simulate for each pattern.

    glider_rules_found = []

    # 1. Iterate through all 256 ECA rules.
    for rule_num in range(256):
        # 2. Filter for compact rules. An ECA is compact if the neighborhood '000'
        # maps to 0. This corresponds to even rule numbers.
        if rule_num % 2 != 0:
            continue

        # Prepare the rule mapping from the rule number.
        # The binary string gives outputs for neighborhoods from '111' down to '000'.
        rule_bin_str = format(rule_num, '08b')
        rule_map = {
            (1, 1, 1): int(rule_bin_str[0]), (1, 1, 0): int(rule_bin_str[1]),
            (1, 0, 1): int(rule_bin_str[2]), (1, 0, 0): int(rule_bin_str[3]),
            (0, 1, 1): int(rule_bin_str[4]), (0, 1, 0): int(rule_bin_str[5]),
            (0, 0, 1): int(rule_bin_str[6]), (0, 0, 0): int(rule_bin_str[7]),
        }

        # 3. For each compact rule, search for a glider.
        found_glider_for_this_rule = False
        # Iterate through possible initial patterns.
        for width in range(1, MAX_INIT_WIDTH + 1):
            # Test all non-trivial patterns of a given width.
            for i in range(1, 2**width):
                initial_pattern_list = [int(b) for b in format(i, f'0{width}b')]
                initial_pattern_tuple = tuple(initial_pattern_list)

                # Setup the simulation grid. It needs to be large enough to
                # contain the pattern's movement without hitting the boundaries.
                grid_size = 2 * MAX_STEPS + width + 4
                # Use two lists as buffers for the current and next grid states
                # to make evolution updates cleaner.
                grid = [0] * grid_size
                
                # Place the pattern in the middle of the grid.
                initial_pos = MAX_STEPS + 2
                grid[initial_pos : initial_pos + width] = initial_pattern_list
                
                # Keep track of the pattern's boundaries for efficiency.
                min_idx, max_idx = initial_pos, initial_pos + width - 1

                # Simulate for MAX_STEPS.
                for step in range(1, MAX_STEPS + 1):
                    next_grid = [0] * grid_size
                    # Evolve the active part of the grid.
                    for j in range(min_idx - 1, max_idx + 2):
                        neighborhood = tuple(grid[j-1 : j+2])
                        next_grid[j] = rule_map[neighborhood]
                    grid = next_grid

                    # Find the boundaries of the new pattern.
                    try:
                        first_one = grid.index(1, min_idx - 2)
                        # Find last '1' by searching the reversed sub-list
                        last_one = (max_idx + 2) - grid[min_idx - 2 : max_idx + 3][::-1].index(1)
                    except ValueError: # Pattern disappeared
                        break

                    # Update boundaries for the next step.
                    min_idx, max_idx = first_one, last_one

                    # Check for boundary collision
                    if min_idx <= 1 or max_idx >= grid_size - 2:
                        break
                        
                    current_pattern_tuple = tuple(grid[min_idx : max_idx + 1])
                    
                    # Check if the pattern is a glider.
                    if current_pattern_tuple == initial_pattern_tuple:
                        displacement = min_idx - initial_pos
                        if displacement != 0:
                            glider_rules_found.append(rule_num)
                            found_glider_for_this_rule = True
                            break
                
                if found_glider_for_this_rule:
                    break
            if found_glider_for_this_rule:
                break
    
    # Print the final result.
    print("Found the following compact ECA rules that have a glider:")
    # Printing the "equation" by showing which numbers were summed.
    equation = " + ".join(map(str, [1]*len(glider_rules_found)))
    print(f"List of rules: {glider_rules_found}")
    print(f"\nCounting them as: {equation} = {len(glider_rules_found)}")
    print("\nTotal number of compact ECAs with a glider:")
    print(len(glider_rules_found))
    
solve()
<<<110>>>