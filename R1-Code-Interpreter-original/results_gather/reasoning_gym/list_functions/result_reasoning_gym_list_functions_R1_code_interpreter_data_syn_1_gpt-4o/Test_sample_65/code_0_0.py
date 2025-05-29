def find_pattern(input_list, output_list):
    differences = []
    for i, o in zip(input_list, output_list):
        differences.append(o - 2 * i)
    return differences

input_examples = [55, 96, 82, 24]
output_examples = [106, 188, 160, 44]

pattern = find_pattern(input_examples, output_examples)
print(pattern)