def rearrange_list(input_list):
    sorted_list = sorted(input_list)
    result = []
    while sorted_list:
        if sorted_list:
            result.append(sorted_list.pop(0))  # smallest
        if sorted_list:
            result.append(sorted_list.pop(-1))  # largest
    return result

input_list = [8, 13, 86, 100, 18, 39, 5, 24]
output_list = rearrange_list(input_list)
print(output_list)