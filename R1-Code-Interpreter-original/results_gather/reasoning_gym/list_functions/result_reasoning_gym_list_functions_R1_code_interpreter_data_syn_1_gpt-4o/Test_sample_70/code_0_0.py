def find_output(input_list):
    smallest = min(input_list)
    if smallest <= 25:
        return [smallest]
    else:
        return [0]

input_list = [11, 7, 89, 24, 78]
output = find_output(input_list)
print(output)