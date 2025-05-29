def analyze_and_predict(input_grid):
    # Get dimensions
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Extract the pattern numbers from the second column onwards
    pattern_numbers = set()
    for i in range(rows):
        for j in range(1, cols):  # Skip first column
            if input_grid[i][j] != input_grid[0][0]:  # Skip the constant number
                pattern_numbers.add(input_grid[i][j])
    
    # Get the smaller and larger numbers
    num1, num2 = min(pattern_numbers), max(pattern_numbers)
    
    # Generate output grid with alternating pattern
    output = []
    for i in range(rows):
        row = []
        start_num = num1 if i % 2 == 0 else num2
        for j in range(cols):
            row.append(start_num if j % 2 == 0 else (num2 if start_num == num1 else num1))
        output.append(row)
    
    # Print the predicted output
    print("Predicted output:")
    for row in output:
        print(" ".join(map(str, row)))

# Test input
test_input = [
    [9, 9, 9, 9, 9, 9],
    [9, 7, 5, 7, 5, 7],
    [9, 5, 7, 5, 7, 5],
    [9, 7, 5, 7, 5, 7],
    [9, 5, 7, 5, 7, 5],
    [9, 7, 5, 7, 5, 7]
]

analyze_and_predict(test_input)