def transform_list(input_list):
    if len(input_list) > 0:
        input_list[0] += 1
        input_list[-1] += 1
    return input_list

input_data = [47, 83, 91, 56, 76, 73, 47]
output_data = transform_list(input_data)
print(output_data)