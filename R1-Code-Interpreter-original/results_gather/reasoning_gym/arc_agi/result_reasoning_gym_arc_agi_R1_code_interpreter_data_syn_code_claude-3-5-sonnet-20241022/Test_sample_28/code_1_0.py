def analyze_pattern(input_grid):
    # Find all non-6 numbers and their positions
    numbers = []
    for i in range(10):
        for j in range(10):
            if input_grid[i][j] != 6:
                numbers.append((i, j, input_grid[i][j]))
    
    # Initialize output grid with 6s
    output = [[6 for _ in range(3)] for _ in range(3)]
    
    # Sort numbers by their position (top to bottom, left to right)
    numbers.sort()
    
    # Map the first 9 numbers to output grid in specific pattern
    # Looking at examples, it seems numbers map in a specific order
    mapping = [(2,0), (2,1), (0,2), (1,1), (2,2), (0,0), (0,1), (1,0), (1,2)]
    
    for idx, (_, _, value) in enumerate(numbers):
        if idx < len(mapping):
            row, col = mapping[idx]
            output[row][col] = value
    
    # Print output grid
    for row in output:
        print(' '.join(map(str, row)))

# Test input
test_input = [
    [6, 6, 0, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 4, 6, 6, 6],
    [3, 6, 6, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 1, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 6, 6, 0, 6],
    [6, 5, 6, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 6, 6, 3, 6],
    [6, 6, 6, 6, 6, 2, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 6, 9, 6, 6]
]

analyze_pattern(test_input)