def process_grid(input_grid):
    height = len(input_grid)
    width = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]
    
    # Process horizontal lines first
    for i in range(height):
        # Find all non-4 numbers in this row
        numbers = {}
        for j in range(width):
            if input_grid[i][j] != 4:
                if input_grid[i][j] not in numbers:
                    numbers[input_grid[i][j]] = []
                numbers[input_grid[i][j]].append(j)
        
        # Fill horizontal lines for each number
        for num, positions in numbers.items():
            if len(positions) >= 2:
                for j in range(min(positions), max(positions) + 1):
                    output_grid[i][j] = num
    
    # Process vertical lines
    for j in range(width):
        # Find all non-4 numbers in this column
        numbers = {}
        for i in range(height):
            if input_grid[i][j] != 4:
                if input_grid[i][j] not in numbers:
                    numbers[input_grid[i][j]] = []
                numbers[input_grid[i][j]].append(i)
        
        # Fill vertical lines for each number
        for num, positions in numbers.items():
            if len(positions) >= 2:
                for i in range(min(positions), max(positions) + 1):
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