def find_output(input_value):
    # Calculate the output as twice the input plus an adjustment
    # The adjustment is a multiple of 6
    base_output = 2 * input_value
    adjustment = (base_output - input_value) // 6 * 6
    output_value = base_output + adjustment
    return output_value

# Test the function with the given input
input_value = 22
output_value = find_output(input_value)
print(output_value)