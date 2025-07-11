def find_ca_rules():
    """
    Identifies the elementary cellular automaton rules that could have produced
    the pattern shown in the image.
    
    The script works by first transcribing the visual pattern into a digital
    grid. It then analyzes the grid to find all unique local transitions that
    occur between generations. These transitions act as constraints. Finally,
    it checks each of the 256 elementary cellular automaton rules against these
    constraints to find the ones that are consistent with the entire observed
    pattern.
    """
    
    # The image shows 8 rows total. This means an initial state (t=0) and
    # 7 subsequent generations (t=1 to t=7). We represent this as 8 rows.
    # The pattern is represented as a grid of 0s (white) and 1s (black).
    # The grid is padded with 0s to handle boundary conditions where the
    # automaton assumes cells outside the visible area are 0.
    # A width of 17 is used, with the initial black cell at index 8.
    grid_str = [
        "00000000100000000", # t=0
        "00000001010000000", # t=1
        "00000010101000000", # t=2
        "00000101010100000", # t=3
        "00001010101010000", # t=4
        "00010101010101000", # t=5
        "00101010101010100", # t=6
        "01010101010101010"  # t=7
    ]
    
    grid = [[int(c) for c in row] for row in grid_str]
    
    # A mapping from a 3-cell neighborhood (left, center, right) to the
    # corresponding bit position in the standard Wolfram rule number.
    # (1,1,1) -> bit 7, (1,1,0) -> bit 6, ..., (0,0,0) -> bit 0
    neighborhood_map = {
        (1, 1, 1): 7, (1, 1, 0): 6, (1, 0, 1): 5, (1, 0, 0): 4,
        (0, 1, 1): 3, (0, 1, 0): 2, (0, 0, 1): 1, (0, 0, 0): 0
    }
    
    # Gather all observed unique transitions (neighborhood -> output) from the grid.
    # This forms a set of constraints that any valid rule must satisfy.
    constraints = set()
    num_rows = len(grid)
    num_cols = len(grid[0])
    
    for t in range(num_rows - 1):
        for x in range(1, num_cols - 1):
            neighborhood = (grid[t][x-1], grid[t][x], grid[t][x+1])
            output = grid[t+1][x]
            constraints.add((neighborhood, output))

    # Test all 256 possible elementary rules.
    possible_rules = []
    for rule_number in range(256):
        is_valid = True
        # Check if the current rule satisfies all observed constraints.
        for neighborhood, expected_output in constraints:
            bit_position = neighborhood_map[neighborhood]
            
            # The rule's output for this neighborhood is the bit at the position.
            rule_output = (rule_number >> bit_position) & 1
            
            if rule_output != expected_output:
                is_valid = False
                break
        
        if is_valid:
            possible_rules.append(rule_number)
            
    # The final answer is the comma-separated list of valid rules, sorted.
    possible_rules.sort()
    print(','.join(map(str, possible_rules)))

find_ca_rules()