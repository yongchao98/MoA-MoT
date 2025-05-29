def find_fives(grid):
    positions = []
    for i in range(len(grid)):
        for j in range(len(grid[i])):
            if grid[i][j] == 5:
                positions.append((i, j))
    return positions

def create_output_grid(input_grid):
    # Create a copy of input grid
    output = [row[:] for row in input_grid]
    
    # Find positions of all 5s
    five_positions = find_fives(input_grid)
    
    # For each 5, place a 4 six steps diagonally
    diagonal_step = 6
    for row, col in five_positions:
        target_row = row + diagonal_step
        target_col = col + diagonal_step
        
        if target_row < len(output) and target_col < len(output[0]):
            output[target_row][target_col] = 4
    
    return output

# Parse input grid
input_str = """5 1 1 1 1 1 1 1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
1 1 1 5 1 1 1 1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
1 1 1 1 1 1 5 1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1"""

# Convert input string to grid
input_grid = [[int(num) for num in row.split()] for row in input_str.split('\n')]

# Generate output grid
output_grid = create_output_grid(input_grid)

# Format output
output_str = '\n'.join(' '.join(str(cell) for cell in row) for row in output_grid)
print(output_str)