def find_numbers_and_positions(grid):
    numbers = []
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if grid[i][j] != 0:
                numbers.append(grid[i][j])
    return sorted(list(set(numbers)))

def create_pattern(number):
    return [
        [0, number, 0],
        [number, number, number],
        [number, 0, number]
    ]

def format_output(patterns):
    result = []
    for pattern in patterns:
        for row in pattern:
            result.append(' '.join(map(str, row)))
    return '\n'.join(result)

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

grid = [[int(x) for x in row.split()] for row in input_str.split('\n')]

# Find all unique non-zero numbers
numbers = find_numbers_and_positions(grid)

# Create patterns for each number
patterns = [create_pattern(num) for num in numbers]

# Print the result
print(format_output(patterns))