import collections

def has_glider(rule_num: int) -> bool:
    """
    Checks if a given Elementary Cellular Automaton rule has a glider.

    It performs a search over small initial configurations and simulates their
    evolution, looking for a pattern that repeats itself at a shifted location.
    The definitions are as follows:
    - An ECA is 'compact' if the rule for '000' is 0 (rule_num is even).
    - A 'glider' is a finite pattern `C` that reappears at a different location
      after some time `t`: Evolve(C, t) = Shift(C, k) for t>0, k!=0.
    """
    # A list where the index is the integer value of the 3-cell neighborhood
    # (e.g., '111'->7, '010'->2) and the value is the next state of the center cell.
    rule_lookup = [(rule_num >> i) & 1 for i in range(8)]

    # Search parameters are chosen to be large enough to find common gliders
    # without taking an excessive amount of time to run.
    MAX_INIT_WIDTH = 8
    MAX_STEPS = 100
    GRID_WIDTH = 2 * MAX_STEPS + MAX_INIT_WIDTH + 20

    # Rule 0 maps everything to 0, so any pattern vanishes. No gliders possible.
    if rule_num == 0:
        return False

    # Iterate through all non-trivial patterns of the given width
    for width in range(1, MAX_INIT_WIDTH + 1):
        for i in range(1, 1 << width):
            init_pattern_list = [int(b) for b in bin(i)[2:].zfill(width)]
            
            grid = [0] * GRID_WIDTH
            start_idx = (GRID_WIDTH - width) // 2
            grid[start_idx : start_idx + width] = init_pattern_list

            # history stores {pattern_shape: (first_time_seen, first_position_seen)}
            history = {}
            
            for t in range(MAX_STEPS + 1):
                # Trim the grid to get the core pattern and its position
                try:
                    first_one = grid.index(1)
                    # Find last '1' by searching the reversed list
                    last_one = len(grid) - 1 - grid[::-1].index(1)
                    current_pattern_tuple = tuple(grid[first_one : last_one + 1])
                    pos = first_one
                except ValueError: # No '1's left in the grid
                    break # The pattern died out

                # Check if this pattern shape has been seen before
                if current_pattern_tuple in history:
                    _prev_t, prev_pos = history[current_pattern_tuple]
                    shift = pos - prev_pos
                    
                    if shift != 0:
                        # Found a glider: the pattern repeated at a new location.
                        return True
                    else:
                        # The pattern is stationary or in a non-moving cycle.
                        # This simulation cannot produce a glider, so we stop it.
                        break
                
                history[current_pattern_tuple] = (t, pos)
                
                # Evolve the grid for the next step
                next_grid = [0] * GRID_WIDTH
                for cell_idx in range(1, GRID_WIDTH - 1):
                    left = grid[cell_idx - 1]
                    center = grid[cell_idx]
                    right = grid[cell_idx + 1]
                    neighborhood_val = (left << 2) | (center << 1) | right
                    next_grid[cell_idx] = rule_lookup[neighborhood_val]
                grid = next_grid

    # If the search completes without finding a glider, we assume none exist
    # within our search limits.
    return False

def solve():
    """
    Finds all compact ECAs that have a glider, prints them, and prints the total count.
    """
    print("Searching for compact ECAs with gliders... (This may take a minute)")
    
    compact_glider_rules = []
    for rule_num in range(256):
        # An ECA is compact iff the rule for '000' is 0, which means rule_num is even.
        if rule_num % 2 == 0:
            if has_glider(rule_num):
                compact_glider_rules.append(rule_num)

    print("\nFound the following compact ECA rules that have a glider:")
    for rule in compact_glider_rules:
        print(rule, end=" ")
    print("\n")
    
    print("The total number of compact ECAs that have a glider is:")
    print(len(compact_glider_rules))

if __name__ == "__main__":
    solve()