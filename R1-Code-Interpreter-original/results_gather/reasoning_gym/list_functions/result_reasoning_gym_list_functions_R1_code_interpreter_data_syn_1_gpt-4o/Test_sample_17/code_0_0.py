def transform_list(input_list):
    if len(input_list) > 1:
        input_list[0] += 1
        input_list[-1] += 1
    return input_list

output = transform_list([93, 18, 13, 55, 94, 49])
print(output)