import collections

def get_grid_from_image():
    """
    Manually transcribes the first few rows of the cellular automaton from the image.
    0 represents a white square, 1 represents a black square.
    The grid is assumed to be 15 cells wide to capture the full pattern.
    """
    grid = [
        [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],  # Row 0 (t=0)
        [0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0],  # Row 1 (t=1)
        [0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0],  # Row 2 (t=2)
    ]
    return grid

def apply_rule(row, rule_map):
    """
    Applies a given ECA rule to a row to generate the next row.
    Assumes fixed '0' (white) boundary conditions.
    """
    width = len(row)
    next_row = [0] * width
    for i in range(width):
        left = row[i - 1] if i > 0 else 0
        center = row[i]
        right = row[i + 1] if i < width - 1 else 0
        neighborhood = (left, center, right)
        next_row[i] = rule_map[neighborhood]
    return next_row

def find_possible_rules():
    """
    Finds all elementary cellular automaton rules that could produce the observed evolution.
    """
    grid = get_grid_from_image()
    initial_row = grid[0]
    target_rows = grid[1:]
    
    possible_rules = []
    
    for rule_number in range(256):
        # Convert rule number to its 8-bit binary representation.
        # This string represents the outputs for neighborhoods 111, 110, ..., 000.
        rule_binary_str = format(rule_number, '08b')
        
        # Create a mapping from neighborhood tuple to output state
        rule_map = {
            (1, 1, 1): int(rule_binary_str[0]),
            (1, 1, 0): int(rule_binary_str[1]),
            (1, 0, 1): int(rule_binary_str[2]),
            (1, 0, 0): int(rule_binary_str[3]),
            (0, 1, 1): int(rule_binary_str[4]),
            (0, 1, 0): int(rule_binary_str[5]),
            (0, 0, 1): int(rule_binary_str[6]),
            (0, 0, 0): int(rule_binary_str[7]),
        }
        
        # Simulate the evolution
        current_row = initial_row
        is_match = True
        for target_row in target_rows:
            next_row = apply_rule(current_row, rule_map)
            if next_row != target_row:
                is_match = False
                break
            current_row = next_row
            
        if is_match:
            possible_rules.append(rule_number)
            
    return sorted(possible_rules)

if __name__ == '__main__':
    solution_rules = find_possible_rules()
    # The problem asks for a comma-separated list
    print(','.join(map(str, solution_rules)))
