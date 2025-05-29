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
    
    # Sort five positions by row number to handle multiple 5s correctly
    five_positions.sort()
    
    # For each 5, place a 4 in the diagonal pattern
    for idx, (row, col) in enumerate(five_positions):
        # Calculate diagonal offset (starts at 2 and increases)
        diagonal_offset = idx + 2
        
        # Check if the target position is within grid bounds
        target_row = row + diagonal_offset
        target_col = col + diagonal_offset
        
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