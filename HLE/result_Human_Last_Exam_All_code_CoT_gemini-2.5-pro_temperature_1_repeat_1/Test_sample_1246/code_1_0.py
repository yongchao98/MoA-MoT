import sys

def find_compact_ecas_with_gliders():
    """
    This script finds all compact Elementary Cellular Automata (ECAs) that possess a glider.

    Method:
    1. An ECA is 'compact' if the '000' neighborhood maps to '0', meaning the rule number is even.
    2. A 'glider' is a finite pattern that repeats itself at a different location after some time steps.
    3. The script iterates through all 128 compact ECAs (even rules from 2 to 254).
    4. For each rule, it tests a set of initial patterns to see if they evolve into a glider.
    5. The search is optimized with parameters sufficient to find known simple gliders:
       - MAX_INITIAL_WIDTH = 8: Tests patterns up to 8 cells wide.
       - MAX_STEPS = 40: Simulates for up to 40 steps.
    6. A rule is counted if any tested initial pattern results in a glider.
    """

    def get_rule_map(rule_num):
        """Converts a rule number (0-255) to a dictionary for fast lookups."""
        rule_map = {}
        # The rule's binary string corresponds to outputs for neighborhoods 111, 110, ..., 000
        binary_rule = format(rule_num, '08b')
        for i in range(8):
            # Neighborhoods from 111 (7) down to 000 (0)
            neighborhood_val = 7 - i
            neighborhood_tuple = tuple(int(x) for x in format(neighborhood_val, '03b'))
            rule_map[neighborhood_tuple] = int(binary_rule[i])
        return rule_map

    def evolve(grid, rule_map):
        """Evolves the grid one step forward according to the rule map."""
        grid_len = len(grid)
        next_grid = [0] * grid_len
        # We assume the infinite grid outside our simulation is all zeros.
        for i in range(grid_len):
            left = grid[i - 1] if i > 0 else 0
            center = grid[i]
            right = grid[i + 1] if i < grid_len - 1 else 0
            neighborhood = (left, center, right)
            next_grid[i] = rule_map.get(neighborhood, 0)
        return next_grid

    def get_pattern_and_pos(grid):
        """Finds the compact pattern (between first/last 1s) and its position."""
        try:
            first_one_idx = grid.index(1)
            # Find last '1' by searching from the end of the list
            last_one_idx = len(grid) - 1 - grid[::-1].index(1)
            # Return pattern as a tuple to be hashable for history set
            pattern = tuple(grid[first_one_idx : last_one_idx + 1])
            return pattern, first_one_idx
        except ValueError:  # No '1's found in grid
            return None, -1

    # --- Simulation Parameters ---
    MAX_INITIAL_WIDTH = 8
    MAX_STEPS = 40
    # Grid must be large enough to contain the pattern's movement and growth
    GRID_WIDTH = MAX_INITIAL_WIDTH + 2 * MAX_STEPS

    # Generate initial patterns to test.
    # We only test "minimal" patterns that start and end with a 1.
    initial_patterns = []
    for width in range(1, MAX_INITIAL_WIDTH + 1):
        if width == 1:
            initial_patterns.append([1])
            continue
        # For a given width, iterate through all possible middle bit combinations
        num_middle_configs = 2**(width - 2)
        for i in range(num_middle_configs):
            middle_str = format(i, f'0{width-2}b')
            pattern_str = '1' + middle_str + '1'
            initial_patterns.append([int(bit) for bit in pattern_str])
    
    glider_rules = []
    
    # Iterate through all compact rules (rule number must be even). Rule 0 is trivial.
    for rule_num in range(2, 256, 2):
        rule_map = get_rule_map(rule_num)
        found_glider_for_this_rule = False

        for initial_pattern_list in initial_patterns:
            initial_pattern_tuple = tuple(initial_pattern_list)
            
            # Set up the simulation grid with the pattern in the center
            grid = [0] * GRID_WIDTH
            start_pos = (GRID_WIDTH - len(initial_pattern_list)) // 2
            grid[start_pos : start_pos + len(initial_pattern_list)] = initial_pattern_list
            
            # History stores (pattern, position) tuples to detect any cycle
            history = set()

            current_grid = grid
            for _ in range(MAX_STEPS):
                current_grid = evolve(current_grid, rule_map)
                pattern, pos = get_pattern_and_pos(current_grid)

                # Condition 1: The pattern died out
                if pattern is None:
                    break
                
                # Condition 2: A glider is found!
                # The pattern is the same as the initial one, but has moved.
                if pattern == initial_pattern_tuple and pos != start_pos:
                    glider_rules.append(rule_num)
                    found_glider_for_this_rule = True
                    break

                # Condition 3: A cycle is detected.
                # If we've seen this exact state (pattern and position) before, stop.
                state = (pattern, pos)
                if state in history:
                    break
                history.add(state)
            
            if found_glider_for_this_rule:
                break  # Move to the next rule

    # --- Final Output ---
    print("Found the following compact ECAs that have a glider:")
    for rule in sorted(glider_rules):
        print(f"Rule {rule}")
    
    print("\nIn total, the number of compact ECAs with a glider is:")
    print(len(glider_rules))

if __name__ == '__main__':
    find_compact_ecas_with_gliders()
<<<20>>>