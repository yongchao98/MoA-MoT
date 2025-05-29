def determine_output(input_list):
    threshold = 90
    for number in input_list:
        if number >= threshold:
            return [1]
    return [0]

input_list = [2, 9, 56, 83, 95, 95]
output = determine_output(input_list)
print(output)