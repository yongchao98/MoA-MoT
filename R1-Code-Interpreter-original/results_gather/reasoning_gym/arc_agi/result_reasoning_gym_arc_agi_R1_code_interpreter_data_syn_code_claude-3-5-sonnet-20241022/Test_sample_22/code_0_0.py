def find_special_position(grid):
    for i in range(len(grid)):
        for j in range(len(grid[i])):
            if grid[i][j] != 1:
                return (i, j), grid[i][j]
    return None, None

def create_output_grid(size, special_pos, special_val):
    output = [[1 for _ in range(size)] for _ in range(size)]
    center_row, center_col = special_pos
    
    # Main diagonal patterns
    # Forward diagonal (top-left to bottom-right)
    r, c = center_row, center_col
    while r >= 0 and c >= 0:
        output[r][c] = special_val
        r -= 1
        c -= 1
    
    r, c = center_row + 1, center_col + 1
    while r < size and c < size:
        output[r][c] = special_val
        r += 1
        c += 1
    
    # Backward diagonal (top-right to bottom-left)
    r, c = center_row, center_col
    while r >= 0 and c < size:
        output[r][c] = special_val
        r -= 1
        c += 1
    
    r, c = center_row + 1, center_col - 1
    while r < size and c >= 0:
        output[r][c] = special_val
        r += 1
        c -= 1
    
    # Additional reflections for larger grids
    if size > 7:
        # Reflect diagonals near the edges
        for i in range(size):
            for j in range(size):
                if output[i][j] == special_val:
                    # Add symmetric reflections
                    if i > size//2 and j < size//2:
                        if i+2 < size and j-2 >= 0:
                            output[i+2][j-2] = special_val
                    if i > size//2 and j > size//2:
                        if i+2 < size and j+2 < size:
                            output[i+2][j+2] = special_val
    
    return output

# Parse input grid
input_grid = [
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 9, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
]

special_pos, special_val = find_special_position(input_grid)
output = create_output_grid(len(input_grid), special_pos, special_val)

# Print the result in the required format
result = ""
for row in output:
    result += " ".join(map(str, row)) + "\n"
print(result.strip())