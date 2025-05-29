def transform_list(input_list):
    output_list = []
    for i, num in enumerate(input_list):
        if i == 0:
            output_list.append(num)
        else:
            output_list.append(num + i)
    return output_list

input_data = [54, 53, 36]
output_data = transform_list(input_data)
print(output_data)