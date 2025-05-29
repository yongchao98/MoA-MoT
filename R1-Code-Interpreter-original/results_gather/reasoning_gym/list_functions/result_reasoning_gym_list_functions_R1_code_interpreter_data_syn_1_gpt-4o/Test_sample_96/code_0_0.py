def transform_list(input_list):
    if len(input_list) > 1:
        input_list[0] += 1
        input_list[-1] += 1
    return input_list

input_data = [67, 46, 61, 42, 64, 42, 68, 45, 53]
output_data = transform_list(input_data)
print(output_data)