def transform_grid(input_grid):
    # Convert input string to list of numbers
    grid = [int(x) for x in input_grid.split()]
    result = grid.copy()
    
    i = 0
    while i < len(grid):
        # Find start of a sequence of 2s
        if grid[i] == 2:
            # Check if we have at least two 2s
            if i + 1 < len(grid) and grid[i + 1] == 2:
                # Transform first two 2s to 8s
                result[i] = 8
                result[i + 1] = 8
                
                # Skip to the end of this sequence of 2s
                i += 2
                while i < len(grid) and grid[i] == 2:
                    i += 1
            else:
                i += 1
        else:
            i += 1
    
    print(' '.join(map(str, result)))

# Test with the given input
test_input = "0 0 2 2 2 2 2 2 2 0 0 0 2 2 2"
transform_grid(test_input)