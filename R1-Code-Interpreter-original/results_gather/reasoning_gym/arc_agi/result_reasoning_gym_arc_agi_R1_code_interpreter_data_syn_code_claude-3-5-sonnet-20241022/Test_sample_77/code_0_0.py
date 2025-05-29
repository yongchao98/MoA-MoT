def print_grid(grid):
    for row in grid:
        print(' '.join(map(str, row)))

def convert_to_grid(input_str):
    rows = []
    numbers = input_str.strip().split()
    for i in range(0, len(numbers), 16):
        rows.append([int(x) for x in numbers[i:i+16]])
    return rows

def process_grid(input_grid):
    n = 16
    output_grid = [row[:] for row in input_grid]  # Create a copy
    
    # Keep first two and last two rows unchanged
    # Process middle section (rows 2 to 13)
    for j in range(16):  # For each column
        for i in range(2, 14):  # For middle rows
            # Check for vertical patterns that need alignment
            if i >= 2 and i <= 5:  # Upper middle section
                if all(input_grid[i][j] == input_grid[i+1][j] for i in range(i, i+2)):
                    output_grid[i][j] = input_grid[i+2][j]
            elif i >= 10 and i <= 13:  # Lower middle section
                if all(input_grid[i][j] == input_grid[i-1][j] for i in range(i, i-2, -1)):
                    output_grid[i][j] = input_grid[i-2][j]
    
    # Mirror the middle vertically
    for i in range(6, 10):
        output_grid[i] = output_grid[15-i]
    
    return output_grid

# Convert input string to grid
input_str = """3 3 1 8 0 1 5 4 4 5 1 0 8 1 3 3
3 8 8 1 1 1 4 4 4 4 1 1 1 8 8 3
1 8 4 0 5 4 4 2 2 2 4 5 0 4 8 1
8 1 0 8 4 4 8 2 2 2 4 4 8 0 1 8
0 1 2 2 2 2 4 2 2 2 3 1 4 5 1 0
1 1 2 2 2 2 4 7 7 4 3 3 4 4 1 1
5 4 2 2 2 2 7 7 7 7 4 4 8 4 4 5
4 4 2 2 2 2 7 3 3 7 7 4 8 8 4 4
4 4 8 8 4 7 7 3 3 7 7 4 8 8 4 4
5 4 4 8 4 4 7 7 7 7 4 4 8 4 4 5
1 1 4 4 3 3 4 7 7 4 3 3 4 4 1 1
0 1 5 4 1 3 4 4 4 4 3 1 4 5 1 0
8 1 0 8 4 4 8 8 8 8 4 4 8 0 1 8
1 8 4 0 5 4 4 8 8 4 4 5 0 4 8 1
3 8 8 1 1 1 4 4 4 4 1 1 1 8 8 3
3 3 1 8 0 1 5 4 4 5 1 0 8 1 3 3"""

input_grid = convert_to_grid(input_str)
output_grid = process_grid(input_grid)
print_grid(output_grid)