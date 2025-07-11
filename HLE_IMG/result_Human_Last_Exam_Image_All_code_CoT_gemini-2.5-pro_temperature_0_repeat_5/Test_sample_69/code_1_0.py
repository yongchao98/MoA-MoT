import collections

def solve_cellular_automaton():
    """
    Finds possible elementary cellular automaton rules for a given pattern.
    The pattern is analyzed for all occurring (neighborhood -> output) transitions.
    A rule is considered possible if for every neighborhood, its output is one
    of the outputs observed for that neighborhood in the pattern.
    """
    # The pattern from the image, 0 for white, 1 for black.
    # The grid is padded with 0s to handle boundary conditions.
    grid_str = [
        "00000000000000000",
        "00000000100000000",
        "00000001110000000",
        "00000011011000000",
        "00000110001100000",
        "00001101010110000",
        "00011011011011000",
        "00110110001101100",
        "01101101010110110",
        "00000000000000000", # Padding row
    ]
    
    grid = [[int(c) for c in row] for row in grid_str]
    
    # This dictionary will store for each neighborhood (key), a set of observed outputs (value)
    observed_outputs = collections.defaultdict(set)
    
    # Analyze the grid to find all observed transitions
    # We iterate up to the second to last row of the pattern
    for t in range(1, len(grid) - 2):
        for i in range(1, len(grid[0]) - 1):
            # The neighborhood is the 3 cells in the row above
            neighborhood = (grid[t][i-1], grid[t][i], grid[t][i+1])
            # The output is the cell in the current row
            output = grid[t+1][i]
            observed_outputs[neighborhood].add(output)

    # Find all rules that are consistent with the observed transitions
    possible_rules = []
    for rule_num in range(256):
        is_possible = True
        # Get the 8-bit binary representation of the rule
        rule_bin = format(rule_num, '08b')
        
        # The rule maps neighborhoods to outputs. The standard order is 111, 110, ..., 000
        neighborhoods = [
            (1,1,1), (1,1,0), (1,0,1), (1,0,0),
            (0,1,1), (0,1,0), (0,0,1), (0,0,0)
        ]
        
        for i in range(8):
            neighborhood = neighborhoods[i]
            rule_output = int(rule_bin[i])
            
            # If a neighborhood was observed in the pattern, the rule's output
            # for it must be one of the observed outputs.
            if neighborhood in observed_outputs:
                if rule_output not in observed_outputs[neighborhood]:
                    is_possible = False
                    break
        
        if is_possible:
            possible_rules.append(rule_num)
            
    print(f"The rules that could have produced this behavior are: {','.join(map(str, sorted(possible_rules)))}")

solve_cellular_automaton()