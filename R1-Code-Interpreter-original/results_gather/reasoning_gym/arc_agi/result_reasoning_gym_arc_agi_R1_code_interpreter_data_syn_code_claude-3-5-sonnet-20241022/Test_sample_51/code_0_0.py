def analyze_transformation(input_grid, output_grid):
    input_rows = len(input_grid)
    input_cols = len(input_grid[0])
    
    # For each input cell, look at corresponding 2x2 block in output
    patterns = {}
    for i in range(input_rows):
        for j in range(input_cols):
            input_val = input_grid[i][j]
            out_i, out_j = i*2, j*2
            block = [
                [output_grid[out_i][out_j], output_grid[out_i][out_j+1]],
                [output_grid[out_i+1][out_j], output_grid[out_i+1][out_j+1]]
            ]
            key = input_val
            if key not in patterns:
                patterns[key] = block
    return patterns

# Example 1
input1 = [[0, 9], [9, 0]]
output1 = [[9, 0, 9, 0], [0, 9, 0, 9], [9, 0, 9, 0], [0, 9, 0, 9]]

# Example 2
input2 = [[9, 9, 1, 9], [9, 1, 1, 9], [1, 9, 1, 9], [9, 9, 9, 1]]
output2 = [
    [1, 1, 9, 1, 1, 1, 9, 1],
    [1, 9, 9, 1, 1, 9, 9, 1],
    [9, 1, 9, 1, 9, 1, 9, 1],
    [1, 1, 1, 9, 1, 1, 1, 9],
    [1, 1, 9, 1, 1, 1, 9, 1],
    [1, 9, 9, 1, 1, 9, 9, 1],
    [9, 1, 9, 1, 9, 1, 9, 1],
    [1, 1, 1, 9, 1, 1, 1, 9]
]

# Example 3
input3 = [[9, 6, 9], [9, 9, 6], [9, 6, 9]]
output3 = [
    [6, 9, 6, 6, 9, 6],
    [6, 6, 9, 6, 6, 9],
    [6, 9, 6, 6, 9, 6],
    [6, 9, 6, 6, 9, 6],
    [6, 6, 9, 6, 6, 9],
    [6, 9, 6, 6, 9, 6]
]

# Analyze patterns from examples
patterns1 = analyze_transformation(input1, output1)
patterns2 = analyze_transformation(input2, output2)
patterns3 = analyze_transformation(input3, output3)

# Test input
test_input = [[9, 9, 3, 3], [9, 3, 9, 9], [9, 9, 9, 3], [3, 9, 3, 9]]

# Generate output based on discovered patterns
def generate_output(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output = [[0] * (cols*2) for _ in range(rows*2)]
    
    for i in range(rows):
        for j in range(cols):
            val = input_grid[i][j]
            # For value 9, output pattern is:
            # [1 1]
            # [1 9]
            if val == 9:
                output[i*2][j*2] = 1
                output[i*2][j*2+1] = 1
                output[i*2+1][j*2] = 1
                output[i*2+1][j*2+1] = 9
            # For value 3, output pattern is:
            # [3 3]
            # [3 1]
            elif val == 3:
                output[i*2][j*2] = 3
                output[i*2][j*2+1] = 3
                output[i*2+1][j*2] = 3
                output[i*2+1][j*2+1] = 1
    
    # Print the output in the required format
    for row in output:
        print(' '.join(map(str, row)))

generate_output(test_input)