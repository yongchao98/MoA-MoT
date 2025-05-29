def find_pattern_shape(grid):
    # Find the first cluster of 1's to determine the pattern
    rows = len(grid)
    cols = len(grid[0])
    pattern = []
    
    # First, find a '1' and capture its pattern
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] == 1:
                # Check 3x3 area around this point
                for di in range(-1, 2):
                    row = []
                    for dj in range(-1, 2):
                        ni, nj = i + di, j + dj
                        if 0 <= ni < rows and 0 <= nj < cols and grid[ni][nj] == 1:
                            row.append(1)
                        else:
                            row.append(0)
                    if any(row):  # Only add rows that have 1's
                        pattern.append(row)
    return pattern

def find_numbers(grid):
    # Find all unique non-zero numbers and their positions
    numbers = set()
    for row in grid:
        for val in row:
            if val != 0:
                numbers.add(val)
    return sorted(list(numbers))

def create_output(pattern, number):
    # Create pattern for a specific number
    output = []
    for row in pattern:
        new_row = []
        for val in row:
            if val == 1:
                new_row.append(str(number))
            else:
                new_row.append('0')
        output.append(new_row)
    return output

# Parse input grid
input_str = """0 0 1 1 0 0 0 0 0 0 0
0 1 1 1 0 0 0 0 0 0 0
0 1 0 0 0 0 5 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 7 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 2 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 4 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0"""

# Convert input to grid
grid = [[int(x) for x in row.split()] for row in input_str.split('\n')]

# Find the pattern from 1's
pattern = find_pattern_shape(grid)

# Find all numbers
numbers = find_numbers(grid)

# Create output for each number
final_output = []
for num in numbers:
    output_pattern = create_output(pattern, num)
    final_output.extend(output_pattern)

# Print the result
for row in final_output:
    print(' '.join(row))