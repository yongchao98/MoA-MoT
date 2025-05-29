def rearrange_list(input_list):
    sorted_list = sorted(input_list)
    output_list = []
    while sorted_list:
        if sorted_list:
            output_list.append(sorted_list.pop())  # largest
        if sorted_list:
            output_list.append(sorted_list.pop(0))  # smallest
    return output_list

input_data = [68, 15, 50]
output_data = rearrange_list(input_data)
print(output_data)