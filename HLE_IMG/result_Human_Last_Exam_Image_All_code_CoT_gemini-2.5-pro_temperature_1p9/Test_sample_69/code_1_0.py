import collections

def solve_cellular_automaton():
    """
    Determines the possible elementary cellular automaton rules that could
    generate the given pattern.
    """
    # 1. Digitize the Pattern from the image.
    # 0 = white, 1 = black
    grid = [
        [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0], # t=0
        [0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0], # t=1
        [0,0,0,0,0,0,1,1,0,1,1,0,0,0,0,0,0], # t=2
        [0,0,0,0,0,1,1,0,0,0,1,1,0,0,0,0,0], # t=3
        [0,0,0,0,1,1,0,1,0,1,0,1,1,0,0,0,0], # t=4
        [0,0,0,1,1,0,1,1,0,1,1,0,1,1,0,0,0], # t=5
        [0,0,1,1,0,1,1,0,0,0,1,1,0,1,1,0,0], # t=6
        [0,1,1,0,1,1,0,1,0,1,0,1,1,0,1,1,0], # t=7
        # Row at t=8 is the last output, not a source of transitions
    ]

    # 2. Extract Transitions
    # A dictionary to store the set of outputs for each observed neighborhood
    transitions = collections.defaultdict(set)

    # Iterate through each step of the evolution (from row t to t+1)
    for t in range(len(grid)):
        prev_row = grid[t]
        # We need the next row to find the output. Stop before the last row.
        if t + 1 >= len(grid):
            # This logic should be on grid t=0..7 and next_row=t=1..8
            # The grid provided for evolution is t=0 to t=8.
            # prev_row is t=0..7, next_row is t=1..8
            # This can't be, grid has 9 rows. prev_row can go up to index 7.
            # next_row up to index 8. The image description has a mistake.
            # Let's fix the grid based on the image description.
            break

    fixed_grid = [
        [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,1,1,0,1,1,0,0,0,0,0,0],
        [0,0,0,0,0,1,1,0,0,0,1,1,0,0,0,0,0],
        [0,0,0,0,1,1,0,1,0,1,0,1,1,0,0,0,0],
        [0,0,0,1,1,0,1,1,0,1,1,0,1,1,0,0,0],
        [0,0,1,1,0,1,1,0,0,0,1,1,0,1,1,0,0],
        [0,1,1,0,1,1,0,1,0,1,0,1,1,0,1,1,0],
        [1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1]
    ]

    for t in range(len(fixed_grid) - 1):
        prev_row = fixed_grid[t]
        next_row = fixed_grid[t+1]
        width = len(prev_row)
        for i in range(width):
            # Get the 3-cell neighborhood assuming 0 (white) boundaries
            left = prev_row[i-1] if i > 0 else 0
            center = prev_row[i]
            right = prev_row[i+1] if i < width-1 else 0
            
            neighborhood = (left, center, right)
            output = next_row[i]
            
            transitions[neighborhood].add(output)

    # 3. Identify Constraints and Ambiguities
    # Standard neighborhood order for rule calculation
    neighborhood_order = [
        (1,1,1), (1,1,0), (1,0,1), (1,0,0),
        (0,1,1), (0,1,0), (0,0,1), (0,0,0)
    ]
    
    rule_template = []
    for n in neighborhood_order:
        if n not in transitions or len(transitions[n]) > 1:
            # Mark as ambiguous if not observed or contradictory
            rule_template.append('?')
        else:
            # Bit is determined if observed and consistent
            rule_template.append(str(list(transitions[n])[0]))
            
    # 4. Construct and Enumerate Rules
    possible_rules = []
    num_ambiguous = rule_template.count('?')
    
    # Iterate through all combinations for the ambiguous bits
    for i in range(2**num_ambiguous):
        temp_rule_list = list(rule_template)
        # Binary string for the current combination of ambiguous bits
        fill_bits = format(i, '0' + str(num_ambiguous) + 'b')
        
        # Substitute '?' with bits from the current combination
        fill_idx = 0
        for j in range(len(temp_rule_list)):
            if temp_rule_list[j] == '?':
                temp_rule_list[j] = fill_bits[fill_idx]
                fill_idx += 1
        
        # Convert binary string to integer rule number
        rule_str = "".join(temp_rule_list)
        rule_num = int(rule_str, 2)
        possible_rules.append(rule_num)

    # 5. Format the Output
    possible_rules.sort()
    print(','.join(map(str, possible_rules)))

solve_cellular_automaton()
>>> 86,94