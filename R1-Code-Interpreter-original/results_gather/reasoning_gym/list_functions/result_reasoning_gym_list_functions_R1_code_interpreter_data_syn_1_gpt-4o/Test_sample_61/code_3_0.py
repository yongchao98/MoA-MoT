def find_output(input_value):
    # Calculate the difference between output and twice the input for known examples
    known_examples = [(46, 88), (16, 28), (10, 16), (76, 148)]
    differences = [output - 2 * input for input, output in known_examples]
    
    # Calculate the average difference
    average_difference = sum(differences) / len(differences)
    
    # Calculate the output using twice the input plus the average difference
    output_value = 2 * input_value + round(average_difference)
    return output_value

# Test the function with the given input
input_value = 22
output_value = find_output(input_value)
print(output_value)