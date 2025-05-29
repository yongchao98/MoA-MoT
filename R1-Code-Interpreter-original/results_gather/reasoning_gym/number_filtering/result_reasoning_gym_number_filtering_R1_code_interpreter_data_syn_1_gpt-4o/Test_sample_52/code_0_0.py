# Original list of numbers as strings
numbers = ['83', '-2.1732', '98.69', '66.622', '-0.2551', '50.86', '21.7205']

# Convert strings to floats and filter out numbers smaller than 93.34
filtered_numbers = [num for num in numbers if float(num) >= 93.34]

# Print the filtered list
print(filtered_numbers)