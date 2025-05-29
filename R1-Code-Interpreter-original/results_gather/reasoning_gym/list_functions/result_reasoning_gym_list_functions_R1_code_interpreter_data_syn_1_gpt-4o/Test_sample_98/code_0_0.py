def find_output(input_list):
    smallest = min(input_list)
    largest = max(input_list)
    difference = largest - smallest
    return difference

input_list = [93, 49, 51, 11, 85, 6, 71, 71, 68, 69]
output = find_output(input_list)
print(output)