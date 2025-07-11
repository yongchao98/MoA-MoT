import collections

def find_compact_ecas_with_gliders():
    """
    Identifies and counts compact Elementary Cellular Automata (ECAs) that possess gliders.

    The script iterates through all 256 ECA rules. It first filters for "compact" rules,
    which are rules that keep finite patterns finite. This property holds only for even-numbered rules.

    For each compact rule, it then searches for a "glider". A glider is a pattern that repeats
    itself after some time but at a new, shifted location. The search is conducted by simulating
    simple starting patterns (a single '1' and small solid blocks of '1's) and checking if their
    evolution leads to a state that has appeared before at a different position.
    """
    
    # Simulation parameters are chosen to be practical for a script to run,
    # while being robust enough to discover many common gliders.
    MAX_STEPS = 512
    GRID_SIZE = 1200
    MAX_INIT_WIDTH = 12

    found_glider_rules = []

    for rule_num in range(256):
        # A rule is "compact" if it maps (0,0,0) to 0. This means the rule number must be even.
        if rule_num % 2 != 0:
            continue

        rule_bits = [(rule_num >> i) & 1 for i in range(8)]
        
        rule_has_glider = False
        
        # We test a set of simple, common initial patterns (seeds).
        for width in range(1, MAX_INIT_WIDTH + 1):
            if rule_has_glider:
                break

            grid = [0] * GRID_SIZE
            start_pos = (GRID_SIZE - width) // 2
            for i in range(width):
                grid[start_pos + i] = 1

            # History stores {pattern_tuple: (time, offset)} to detect repeats.
            history = {}

            for t in range(MAX_STEPS):
                # Extract the compact pattern and its offset from the grid.
                try:
                    first_one = grid.index(1)
                    # Find the last '1' by searching from the end of the list.
                    last_one = GRID_SIZE - 1 - grid[::-1].index(1)
                    current_pattern_tuple = tuple(grid[first_one:last_one + 1])
                    current_offset = first_one
                except ValueError:
                    # The configuration has no more '1's and has died out.
                    break
                
                # If this exact pattern has been seen before...
                if current_pattern_tuple in history:
                    _prev_time, prev_offset = history[current_pattern_tuple]
                    # ...check if its location is different.
                    if current_offset != prev_offset:
                        rule_has_glider = True
                        break
                
                history[current_pattern_tuple] = (t, current_offset)
                
                # Evolve the grid for one time step.
                next_grid = [0] * GRID_SIZE
                for i in range(1, GRID_SIZE - 1):
                    neighborhood_idx = (grid[i - 1] << 2) | (grid[i] << 1) | grid[i + 1]
                    next_grid[i] = rule_bits[neighborhood_idx]
                grid = next_grid

                # Stop if the pattern grows to the edge of our simulation grid.
                if grid[0] == 1 or grid[GRID_SIZE - 1] == 1:
                    break
        
        if rule_has_glider:
            found_glider_rules.append(rule_num)

    print(f"Based on the simulation, the number of compact ECAs with a glider is {len(found_glider_rules)}.")
    print("The identified rules are:")
    print(found_glider_rules)
    print(f"\n<<< {len(found_glider_rules)} >>>")


find_compact_ecas_with_gliders()