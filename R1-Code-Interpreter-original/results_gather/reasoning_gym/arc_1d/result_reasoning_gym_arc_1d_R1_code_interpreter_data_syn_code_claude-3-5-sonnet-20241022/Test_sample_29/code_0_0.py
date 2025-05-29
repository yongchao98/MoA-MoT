def analyze_and_transform(input_grid):
    # Convert input string to list of numbers
    grid = [int(x) for x in input_grid.split()]
    result = grid.copy()
    
    # Scan for consecutive pairs of 2s
    i = 0
    while i < len(grid) - 1:
        if grid[i] == 2 and grid[i+1] == 2:
            # Check if this is part of a larger sequence of 2s
            if (i == 0 or grid[i-1] != 2) and (i+2 == len(grid) or grid[i+2] != 2):
                # If it's exactly a pair of 2s, transform to 8s
                result[i] = 8
                result[i+1] = 8
            i += 2
        else:
            i += 1
    
    # Print result
    print(' '.join(map(str, result)))

# Test with the given input
test_input = "0 0 2 2 2 2 2 2 2 0 0 0 2 2 2"
analyze_and_transform(test_input)