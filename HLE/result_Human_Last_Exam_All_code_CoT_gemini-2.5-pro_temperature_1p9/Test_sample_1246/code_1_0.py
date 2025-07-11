import numpy as np
import sys
from contextlib import redirect_stdout
import io

def has_glider_for_rule(rule_num, max_width, max_steps, grid_size):
    """
    Checks if a given ECA rule has a glider by searching through small initial patterns.

    Args:
        rule_num (int): The ECA rule number (0-255).
        max_width (int): The maximum width of initial patterns to test.
        max_steps (int): The maximum number of simulation steps.
        grid_size (int): The size of the simulation grid.

    Returns:
        bool: True if a glider is found, False otherwise.
    """
    # A compact ECA must map the '000' neighborhood to 0. This means the rule number must be even.
    if rule_num % 2 != 0:
        return False
        
    # Rule 0 maps everything to 0, so no non-trivial configuration can persist.
    if rule_num == 0:
        return False

    # Create a vector for the rule map for efficient numpy lookup.
    # The neighborhood (l, c, r) is treated as a binary number lcr, which gives an index from 0 to 7.
    # The rule number's binary representation b7b6...b0 corresponds to outputs for neighborhoods 111, 110, ..., 000.
    # So we need to map index i to bit b_i, which requires reversing the binary string.
    rule_bin = np.binary_repr(rule_num, 8)
    rule_map_vector = np.array([int(bit) for bit in rule_bin[::-1]], dtype=np.int8)

    checked_initial_cores = set()

    for width in range(1, max_width + 1):
        # Iterate through all 2^width - 1 non-trivial patterns of a given width
        for i in range(1, 1 << width):
            initial_pattern_list = [int(bit) for bit in np.binary_repr(i, width)]

            # Determine the core pattern and avoid re-checking symmetrical or padded versions
            non_zero_indices = [j for j, val in enumerate(initial_pattern_list) if val == 1]
            first_one, last_one = non_zero_indices[0], non_zero_indices[-1]
            initial_core_tuple = tuple(initial_pattern_list[first_one : last_one + 1])
            
            if initial_core_tuple in checked_initial_cores:
                continue
            checked_initial_cores.add(initial_core_tuple)

            # Set up the simulation grid with the pattern in the center
            grid = np.zeros(grid_size, dtype=np.int8)
            core_width = len(initial_core_tuple)
            start_pos = (grid_size - core_width) // 2
            grid[start_pos : start_pos + core_width] = initial_core_tuple
            
            initial_core_pos = start_pos

            for _ in range(max_steps):
                # Evolve one step
                # We pad the grid to handle boundaries as quiescent (all zeros)
                padded_grid = np.pad(grid, 1, 'constant')
                # Calculate neighborhood indices: left*4 + center*2 + right*1
                indices = (padded_grid[:-2] * 4 + padded_grid[1:-1] * 2 + padded_grid[2:])
                grid = rule_map_vector[indices]

                # Find the core of the current pattern
                current_non_zero = np.where(grid == 1)[0]
                if len(current_non_zero) == 0: # Pattern died out
                    break 
                
                first_one, last_one = current_non_zero[0], current_non_zero[-1]
                current_core_tuple = tuple(grid[first_one : last_one + 1])
                
                # Check for glider condition: same core pattern, different location
                if current_core_tuple == initial_core_tuple:
                    shift = first_one - initial_core_pos
                    if shift != 0:
                        return True
                        
    return False

def find_compact_ecas_with_gliders():
    """
    Main function to find and count all compact ECAs with gliders.
    """
    # Parameters for the search. These are chosen to be practical for a script.
    # A deeper search might yield more results but take significantly longer.
    MAX_INITIAL_WIDTH = 10
    MAX_STEPS = 256
    GRID_SIZE = 256

    glider_rules = []
    
    # Only even-numbered rules can be compact
    for rule_num in range(0, 256, 2):
        # Redirect stdout to hide potential progress indicators within the function
        f = io.StringIO()
        with redirect_stdout(f):
             # To show progress during the long search, print dots to stderr
            print(f"Checking Rule {rule_num}...", file=sys.stderr)
            if has_glider_for_rule(rule_num, MAX_INITIAL_WIDTH, MAX_STEPS, GRID_SIZE):
                glider_rules.append(rule_num)
    
    print("Based on the computational search, the compact ECAs found to have gliders are:")
    # To format the "equation" as requested
    equation_str = " + ".join(["1"] * len(glider_rules))
    print(f"List of rules: {glider_rules}")
    # print(f"This gives the equation: {equation_str} = {len(glider_rules)}")
    print(f"\nTotal number of compact ECAs with a glider: {len(glider_rules)}")
    return len(glider_rules)

if __name__ == '__main__':
    count = find_compact_ecas_with_gliders()
    # The final answer is enclosed in <<<>>>
    print(f"\n<<<{count}>>>")
