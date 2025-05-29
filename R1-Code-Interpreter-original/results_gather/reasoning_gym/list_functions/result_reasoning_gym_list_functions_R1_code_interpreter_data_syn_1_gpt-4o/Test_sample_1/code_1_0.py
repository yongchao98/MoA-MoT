input_list = [36, 41, 72, 98, 88, 32, 16, 8, 41]
output_list = [num for num in input_list if num < 50]
print([min(output_list)] if output_list else [0])