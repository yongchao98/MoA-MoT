def transform_list(input_list):
    if len(input_list) > 1:
        input_list[0] += 1
        input_list[-1] += 1
    return input_list

input_data = [81, 39, 53, 6, 60]
output_data = transform_list(input_data)
print(output_data)