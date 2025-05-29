def transform_grid(input_grid):
    # Create output grid as copy of input
    output = [row[:] for row in input_grid]
    
    # Find the two special numbers
    special_numbers = []
    for j in range(10):
        if input_grid[1][j] != 3:
            special_numbers.append(input_grid[1][j])
        if input_grid[8][j] != 3:
            special_numbers.append(input_grid[8][j])
    
    # If we don't have exactly 2 special numbers, something's wrong
    if len(special_numbers) != 2:
        return None
    
    num1, num2 = special_numbers[0], special_numbers[1]
    
    # Fill the inner part of the 9x9 border
    # Top-left quadrant
    for i in range(3, 5):
        for j in range(3, 5):
            output[i][j] = num1
    
    # Top-right quadrant
    for i in range(3, 5):
        for j in range(5, 7):
            output[i][j] = num2
    
    # Bottom-left quadrant
    for i in range(5, 7):
        for j in range(3, 5):
            output[i][j] = num2
    
    # Bottom-right quadrant
    for i in range(5, 7):
        for j in range(5, 7):
            output[i][j] = num1
    
    # Convert to string format
    result = ""
    for row in output:
        result += " ".join(map(str, row)) + "\n"
    
    print(result.strip())

# Test input
test_input = [
    [3,3,3,3,3,3,3,3,3,3],
    [3,8,3,3,3,3,3,3,1,3],
    [3,3,9,9,9,9,9,9,3,3],
    [3,3,9,3,3,3,3,9,3,3],
    [3,3,9,3,3,3,3,9,3,3],
    [3,3,9,3,3,3,3,9,3,3],
    [3,3,9,3,3,3,3,9,3,3],
    [3,3,9,9,9,9,9,9,3,3],
    [3,3,3,3,3,3,3,3,3,3],
    [3,3,3,3,3,3,3,3,3,3]
]

transform_grid(test_input)