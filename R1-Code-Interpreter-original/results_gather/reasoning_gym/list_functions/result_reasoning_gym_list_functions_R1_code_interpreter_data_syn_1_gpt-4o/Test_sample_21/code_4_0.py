def find_transformed_output(input_list):
    total_sum = sum(input_list)
    # Try different transformations to match the pattern
    transformed_output = total_sum % 100  # Example transformation
    return transformed_output

input_list = [66, 85, 52, 8, 43, 36, 20, 12]
output = find_transformed_output(input_list)
print(output)