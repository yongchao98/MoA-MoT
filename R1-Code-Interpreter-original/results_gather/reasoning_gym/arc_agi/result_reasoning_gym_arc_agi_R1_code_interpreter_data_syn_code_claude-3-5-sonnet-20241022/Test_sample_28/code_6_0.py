def find_output_pattern(input_grid):
    # Find all non-6 numbers and their positions
    numbers = []
    for i in range(10):
        for j in range(10):
            if input_grid[i][j] != 6:
                numbers.append((i, j, input_grid[i][j]))
    
    # Sort numbers by column position (right to left)
    numbers.sort(key=lambda x: (-x[1]))
    
    # Initialize output grid
    output = [[0 for _ in range(3)] for _ in range(3)]
    
    # First three numbers (from rightmost columns) go to right column, top to bottom
    right_nums = numbers[:3]
    right_nums.sort(key=lambda x: x[0])  # sort by row
    for i, (_, _, val) in enumerate(right_nums):
        output[i][2] = val
        
    # Next three numbers go to middle column
    mid_nums = numbers[3:6]
    mid_nums.sort(key=lambda x: x[0])
    for i, (_, _, val) in enumerate(mid_nums):
        output[i][1] = val
        
    # Remaining numbers go to left column
    left_nums = numbers[6:]
    left_nums.sort(key=lambda x: x[0])
    for i, (_, _, val) in enumerate(left_nums):
        output[i][0] = val
    
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

find_output_pattern(test_input)