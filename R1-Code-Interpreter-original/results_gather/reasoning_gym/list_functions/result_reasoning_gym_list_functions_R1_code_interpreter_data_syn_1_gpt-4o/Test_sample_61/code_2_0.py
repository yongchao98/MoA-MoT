def find_output(input_value):
    # Calculate the ratio of output to input for known examples
    known_examples = [(46, 88), (16, 28), (10, 16), (76, 148)]
    ratios = [output / input for input, output in known_examples]
    
    # Calculate the average ratio
    average_ratio = sum(ratios) / len(ratios)
    
    # Calculate the output using the average ratio
    output_value = round(input_value * average_ratio)
    return output_value

# Test the function with the given input
input_value = 22
output_value = find_output(input_value)
print(output_value)