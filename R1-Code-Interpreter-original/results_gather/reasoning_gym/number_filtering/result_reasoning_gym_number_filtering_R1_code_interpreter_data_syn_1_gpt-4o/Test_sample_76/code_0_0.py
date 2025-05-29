# Original list of numbers as strings
numbers = ['43', '76.088', '45.655', '-61', '-40.1']

# Convert strings to floats and filter numbers larger than -33.31
filtered_numbers = [num for num in numbers if float(num) > -33.31]

# Print the filtered list
print(filtered_numbers)