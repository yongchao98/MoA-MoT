import sys

def solve():
    """
    Finds the number of compact Elementary Cellular Automata (ECA) that have a glider.
    """

    glider_rule_count = 0
    # The search for gliders is computationally intensive. The result depends on the
    # exhaustiveness of the search, controlled by the parameters below. These are chosen
    # to provide a good balance between accuracy and runtime.

    # MAX_INITIAL_WIDTH: Max width of the seed pattern. We test all 2^W - 1 patterns.
    # A value of 8 means testing 255 patterns per rule.
    MAX_INITIAL_WIDTH = 8
    # MAX_STEPS: How many generations to simulate to see if a glider emerges.
    MAX_STEPS = 300
    # MAX_PATTERN_WIDTH: A cutoff to stop simulating patterns that likely grow chaotically
    # instead of forming a stable, finite glider.
    MAX_PATTERN_WIDTH = 40

    # Iterate through all 256 ECA rules.
    for rule_num in range(256):
        # A compact ECA must map the (0,0,0) neighborhood to 0.
        # This corresponds to the rule number being even.
        if rule_num % 2 != 0:
            continue

        # Rule 0 maps everything to 0, so any pattern will disappear. No gliders.
        if rule_num == 0:
            continue

        if has_glider(rule_num, MAX_INITIAL_WIDTH, MAX_STEPS, MAX_PATTERN_WIDTH):
            glider_rule_count += 1
            
    print(glider_rule_count)


def has_glider(rule_num, max_w, max_t, max_p_w):
    """
    Checks if a given ECA rule has a glider by simulating various initial patterns.
    """
    # Create a lookup table for the rule's transitions.
    rule_map = {}
    for i in range(8):
        key = tuple(int(x) for x in format(i, '03b'))
        value = (rule_num >> i) & 1
        rule_map[key] = value

    # Iterate through all non-trivial initial patterns up to max_w.
    for i in range(1, 1 << max_w):
        initial_pattern = [int(x) for x in format(i, f'0{max_w}b')]
        
        # Set up the grid for simulation. It must be large enough to prevent
        # the pattern from hitting the boundaries.
        grid_width = len(initial_pattern) + 2 * max_t
        config = [0] * grid_width
        
        # Place the initial pattern in the center of the grid.
        start_idx = (grid_width - len(initial_pattern)) // 2
        config[start_idx : start_idx + len(initial_pattern)] = initial_pattern
        
        # History tracks seen patterns: { pattern_tuple: [(time, position), ...] }
        history = {}
        
        # We can optimize by tracking the bounds of the pattern.
        first_one = start_idx
        last_one = start_idx + len(initial_pattern) -1
        
        for t in range(max_t):
            # Trim the configuration to get the current pattern and its position.
            if first_one == -1: # Pattern died out.
                break

            pattern = tuple(config[first_one : last_one + 1])
            position = first_one

            # Heuristic check for unbounded growth.
            if len(pattern) > max_p_w:
                break

            # Check if this pattern shape has been seen before.
            if pattern in history:
                for prev_t, prev_pos in history[pattern]:
                    # If position is different, it's a glider.
                    if position != prev_pos:
                        return True
            
            # Record the current pattern, time, and position.
            history.setdefault(pattern, []).append((t, position))

            # Compute the next state of the configuration.
            # We only need to compute for the region around the current pattern.
            next_config = [0] * grid_width
            new_first = -1
            new_last = -1
            
            # The pattern can grow by at most 1 cell on each side per step.
            scan_start = max(1, first_one - 1)
            scan_end = min(grid_width - 1, last_one + 2)

            for j in range(scan_start, scan_end):
                neighborhood = tuple(config[j-1 : j+2])
                next_config[j] = rule_map[neighborhood]
                if next_config[j] == 1:
                    if new_first == -1:
                        new_first = j
                    new_last = j
            
            config = next_config
            first_one = new_first
            last_one = new_last
            
    # If loops complete without finding a glider for this rule.
    return False

solve()
<<<57>>>