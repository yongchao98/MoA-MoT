def transform_list(input_list):
    if not input_list:
        return input_list
    # Increment the first and last elements by 1
    input_list[0] += 1
    input_list[-1] += 1
    return input_list

# Given input
input_data = [68, 43, 92, 98, 12, 51, 91]
# Transform the input data
output_data = transform_list(input_data)
print(output_data)