def calculate_difference_sum(input_list):
    difference_sum = 0
    for i in range(len(input_list) - 1):
        difference_sum += abs(input_list[i] - input_list[i + 1])
    return difference_sum

input_list = [66, 85, 52, 8, 43, 36, 20, 12]
output = calculate_difference_sum(input_list)
print(output)