def find_pattern(input_list, output_list):
    # Check if the output is a sublist of the input
    for start in range(len(input_list)):
        for end in range(start + 1, len(input_list) + 1):
            if input_list[start:end] == output_list:
                return start, end - start
    return None

# Given examples
examples = [
    ([75, 4, 30, 4, 50, 47, 9, 55], [75, 4, 30, 4]),
    ([44, 13, 75, 66, 20, 13, 36], [75, 66, 13]),
    ([79, 1, 100, 20, 97, 68, 12, 8], [1, 100, 97, 68, 12]),
    ([33, 95, 71, 91, 8, 50, 22, 76], [91, 8, 50, 22])
]

# Analyze each example
patterns = [find_pattern(input_list, output_list) for input_list, output_list in examples]
print(patterns)