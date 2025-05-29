def process_grid(input_grid):
    height = len(input_grid)
    width = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]
    
    # Process each unique number (except 4) separately
    unique_numbers = set(num for row in input_grid for num in row if num != 4)
    
    for num in unique_numbers:
        # Process horizontal lines
        for i in range(height):
            positions = [j for j in range(width) if input_grid[i][j] == num]
            if len(positions) >= 2:
                for j in range(positions[0], positions[-1] + 1):
                    output_grid[i][j] = num
        
        # Process vertical lines
        for j in range(width):
            positions = [i for i in range(height) if input_grid[i][j] == num]
            if len(positions) >= 2:
                for i in range(positions[0], positions[-1] + 1):
                    output_grid[i][j] = num
    
    return output_grid

# Test input
test_input = [
    [4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4],
    [4,4,4,7,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4],
    [4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4],
    [4,1,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,1,4],
    [4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4],
    [4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4],
    [4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4],
    [4,4,4,4,4,4,4,8,4,4,4,4,4,4,4,4,4,4,4,4],
    [4,4,4,4,4,4,4,4,4,4,4,4,4,4,0,4,4,4,4,4],
    [4,4,9,4,4,4,4,4,4,9,4,4,4,4,4,4,4,4,4,4],
    [4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4],
    [4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4],
    [4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4],
    [4,4,4,4,4,4,4,8,4,4,4,4,4,4,4,4,4,4,4,4],
    [4,4,4,4,4,4,4,4,4,4,4,4,4,4,0,4,4,4,4,4],
    [4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4],
    [4,4,4,7,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4],
    [4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4],
    [4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4],
    [4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4]
]

# Process the grid and format the output
result = process_grid(test_input)
output = ""
for row in result:
    output += " ".join(map(str, row)) + "\n"
print(output.strip())