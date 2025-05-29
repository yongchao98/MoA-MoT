def transform_grid(input_grid):
    output = input_grid.copy()
    n = len(input_grid)
    
    # Find groups of non-zero numbers and mark adjacent zeros
    i = 0
    while i < n:
        # If we find a non-zero number
        if input_grid[i] != 0:
            # Find the end of this group
            start = i
            while i < n and input_grid[i] != 0:
                i += 1
            end = i - 1
            
            # If the group length is >= 2, mark adjacent zeros with 3
            if end > start:
                # Mark left side
                if start > 0 and input_grid[start-1] == 0:
                    output[start-1] = 3
                # Mark right side
                if end < n-1 and input_grid[end+1] == 0:
                    output[end+1] = 3
        else:
            i += 1
    
    return output

# Test input
test_input = [4, 0, 0, 5, 5, 0, 0, 0, 9, 9, 9, 0, 0, 0, 0]
result = transform_grid(test_input)
print(' '.join(map(str, result)))