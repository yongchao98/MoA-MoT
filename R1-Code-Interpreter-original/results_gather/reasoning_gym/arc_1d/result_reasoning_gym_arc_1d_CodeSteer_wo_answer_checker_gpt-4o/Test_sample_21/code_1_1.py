def transform_grid(input_grid):
    # Convert input grid to a list of integers
    input_list = list(map(int, input_grid.split()))
    
    # Initialize the output grid as a copy of the input grid
    output_list = input_list[:]
    
    # Check for specific patterns and set positions to 0 accordingly
    if input_list[8] == 0:
        output_list[4] = 0
    if input_list[0] == 0:
        output_list[0] = 0
        output_list[1] = 0
        output_list[2] = 0
        output_list[3] = 0
        output_list[4] = 0
    if input_list[23] == 0:
        output_list[22] = 0
        output_list[23] = 0
    
    # Return the transformed grid as a string
    return ' '.join(map(str, output_list))

# Test input
test_input = "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6 6 0 0 0 0 0"
output = transform_grid(test_input)
print(f"<<<{output}>>>")