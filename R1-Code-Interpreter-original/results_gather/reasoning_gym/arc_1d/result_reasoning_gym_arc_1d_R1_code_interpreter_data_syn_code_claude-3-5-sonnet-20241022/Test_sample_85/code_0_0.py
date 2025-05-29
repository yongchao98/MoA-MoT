def find_groups(grid):
    groups = []
    current_group = []
    
    for i, num in enumerate(grid):
        if num != 0:
            current_group.append((i, num))
        elif current_group:
            groups.append(current_group)
            current_group = []
    
    if current_group:
        groups.append(current_group)
    
    return groups

def apply_transformation(input_grid):
    # Convert input string to list of integers
    grid = [int(x) for x in input_grid.split()]
    
    # Find groups of non-zero numbers
    groups = find_groups(grid)
    
    # Create output grid filled with zeros
    output = [0] * len(grid)
    
    # For each group, move it to the left
    for group in groups:
        # Get the numbers from the group
        numbers = [x[1] for x in group]
        # Place them at the leftmost available position
        start_pos = 0
        while start_pos < len(output) and output[start_pos] != 0:
            start_pos += 1
        
        for i, num in enumerate(numbers):
            if start_pos + i < len(output):
                output[start_pos + i] = num
    
    return ' '.join(map(str, output))

# Test input
test_input = "0 0 0 0 0 0 0 0 0 0 0 0 7 0 0"
print(apply_transformation(test_input))