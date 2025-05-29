def find_specific_sixes(grid):
    first_six = None
    second_half_six = None
    last_six = None
    rows = len(grid)
    cols = len(grid[0])
    
    # Find the first occurrence of 6
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] == 6:
                if first_six is None:
                    first_six = (i, j)
                # Check if it's in the second half of the grid
                if i >= rows // 2 and second_half_six is None:
                    second_half_six = (i, j)
                # Always update the last occurrence
                last_six = (i, j)
    
    return [first_six, second_half_six, last_six]

# Test input grid
input_grid = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0],
    [0, 3, 0, 0, 0, 0, 0, 6, 0, 0, 0, 3, 0, 0, 0],
    [0, 3, 0, 0, 0, 3, 0, 6, 0, 0, 0, 3, 0, 0, 0],
    [0, 3, 0, 0, 0, 3, 0, 6, 0, 0, 0, 3, 0, 0, 0],
    [0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 6, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 6, 0],
    [0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0, 3, 0, 6, 0],
    [0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
]

output_grid = find_specific_sixes(input_grid)
print(output_grid)