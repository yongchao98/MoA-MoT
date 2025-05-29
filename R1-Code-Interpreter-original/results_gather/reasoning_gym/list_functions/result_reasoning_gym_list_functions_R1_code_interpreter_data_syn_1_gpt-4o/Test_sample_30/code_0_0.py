def transform_list(input_list):
    output_list = [input_list[0]]  # First element remains unchanged
    for i in range(1, len(input_list)):
        output_list.append(input_list[i] + i)
    return output_list

input_data = [89, 28, 92, 11, 42]
output_data = transform_list(input_data)
print(output_data)