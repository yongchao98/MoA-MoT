def rearrange_list(input_list):
    result = []
    input_list.sort()
    while input_list:
        result.append(input_list.pop(0))  # Smallest
        if input_list:
            result.append(input_list.pop(-1))  # Largest
    return result

input_list = [24, 63, 19, 7, 36, 91, 78, 81, 69, 55]
output_list = rearrange_list(input_list)
print(output_list)