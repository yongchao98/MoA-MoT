# Original list of numbers as strings
numbers = ['-63.7', '61.95', '-41', '-2.77', '-86.975', '-34.0', '-33', '79', '26.96', '-94.7']

# Convert strings to floats and filter out numbers smaller than -72.485
filtered_numbers = [num for num in numbers if float(num) >= -72.485]

# Print the filtered list
print(filtered_numbers)