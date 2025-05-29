def process_grid(input_grid):
    height = len(input_grid)
    width = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]  # Create a copy of input grid
    
    # Process each row
    for i in range(height):
        # Find all non-4 numbers in the row and their positions
        numbers = {}
        for j in range(width):
            if input_grid[i][j] != 4:
                if input_grid[i][j] not in numbers:
                    numbers[input_grid[i][j]] = []
                numbers[input_grid[i][j]].append(j)
        
        # For each unique number, check if there are multiple instances
        for num, positions in numbers.items():
            # Also check vertically for the same number
            for pos in positions:
                vertical_positions = []
                for row in range(height):
                    if input_grid[row][pos] == num:
                        vertical_positions.append(row)
                
                # If we have multiple instances vertically
                if len(vertical_positions) > 1:
                    min_pos = min(vertical_positions)
                    max_pos = max(vertical_positions)
                    # Fill vertically
                    for row in range(min_pos, max_pos + 1):
                        output_grid[row][pos] = num
                
            # If we have multiple instances horizontally
            if len(positions) > 1:
                min_pos = min(positions)
                max_pos = max(positions)
                # Fill horizontally
                for j in range(min_pos, max_pos + 1):
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