import time

def solve():
    """
    This script finds all Elementary Cellular Automata (ECA) that are 'compact' and have a 'glider'.

    The plan is as follows:
    1.  Identify all 'compact' ECAs. An ECA is compact if the rule for the '000' neighborhood is 0,
        which means the rule number must be even. There are 128 such rules.
    2.  For each of these 128 rules, search for a 'glider'. A glider is a finite pattern of 1s that
        reappears in a different location after some time.
    3.  The search is done by simulation. We test various small initial patterns (up to a certain width)
        and evolve them for a set number of steps.
    4.  If, during the simulation, a pattern reappears identical to its start state but shifted,
        we've found a glider for that rule. We then count that rule and move to the next.
    5.  The final count is the answer to the problem.
    """
    # Simulation parameters are chosen to be reasonably thorough yet run in a feasible amount of time.
    W_MAX = 8      # Maximum width of initial patterns to test
    T_MAX = 300    # Maximum simulation steps for each pattern
    GRID_SIZE = 600  # Size of the simulation grid (must be large enough to avoid boundary effects)

    def has_glider_for_rule(rule):
        """Checks a single ECA rule for a glider by simulating small initial patterns."""
        rule_bits = [(rule >> i) & 1 for i in range(8)]

        # Generate a list of non-trivial initial patterns to test
        initial_patterns = []
        for width in range(1, W_MAX + 1):
            if width == 1:
                initial_patterns.append((1,))
                continue
            # Patterns must start and end with 1, so we vary the middle w-2 bits
            num_middle_patterns = 1 << (width - 2)
            for i in range(num_middle_patterns):
                middle = tuple(int(b) for b in bin(i)[2:].zfill(width - 2))
                pattern = (1,) + middle + (1,)
                initial_patterns.append(pattern)

        # Test each pattern for the given rule
        for initial_pattern_tuple in initial_patterns:
            grid = [0] * GRID_SIZE
            pattern_len = len(initial_pattern_tuple)
            start_pos = (GRID_SIZE - pattern_len) // 2
            
            grid[start_pos : start_pos + pattern_len] = initial_pattern_tuple
            
            # Run the simulation for T_MAX steps
            for _ in range(1, T_MAX + 1):
                next_grid = [0] * GRID_SIZE
                # Apply the ECA rule to each cell
                for i in range(1, GRID_SIZE - 1):
                    neighborhood_value = (grid[i-1] << 2) | (grid[i] << 1) | grid[i+1]
                    next_grid[i] = rule_bits[neighborhood_value]
                grid = next_grid

                # Find the bounding box of the current pattern
                try:
                    left_pos = grid.index(1)
                    right_pos = GRID_SIZE - 1 - grid[::-1].index(1)
                except ValueError:
                    # The pattern has died out (all zeros)
                    break  

                # If pattern is too close to the boundary, simulation is unreliable
                if left_pos < 5 or right_pos > GRID_SIZE - 6:
                    break

                # Extract the current pattern shape
                current_pattern_tuple = tuple(grid[left_pos : right_pos + 1])
                
                # Check if the shape matches the initial one
                if current_pattern_tuple == initial_pattern_tuple:
                    displacement = left_pos - start_pos
                    # If it has moved, it's a glider
                    if displacement != 0:
                        return True
                        
        return False # No glider found for this rule with the tested patterns

    glider_rules = []
    # Iterate through all 128 compact ECA rules (even numbers from 0 to 254)
    for rule in range(0, 256, 2):
        if has_glider_for_rule(rule):
            glider_rules.append(rule)

    print(f"Searching all 128 compact ECAs for gliders...")
    print(f"Found {len(glider_rules)} compact ECAs that have a glider.")
    print("These are the rule numbers:")
    print(sorted(glider_rules))
    print(f"<<<{len(glider_rules)}>>>")

if __name__ == '__main__':
    solve()