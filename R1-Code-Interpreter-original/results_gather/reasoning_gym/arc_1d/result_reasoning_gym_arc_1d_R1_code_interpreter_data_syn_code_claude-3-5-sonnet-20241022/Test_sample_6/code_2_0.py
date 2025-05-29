def generate_output(input_grid):
    base_num = input_grid[0]
    output = [0] * 23
    
    # First row: positions 0-4
    for i in range(5):
        output[i] = base_num
    
    # Second row: positions 6-10
    for i in range(6, 11):
        output[i] = base_num
    
    # Third row: positions 13-17
    for i in range(13, 18):
        output[i] = base_num
    
    return output

# Test with the given input
test_input = [1,1,1,1,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
result = generate_output(test_input)
print(' '.join(map(str, result)))