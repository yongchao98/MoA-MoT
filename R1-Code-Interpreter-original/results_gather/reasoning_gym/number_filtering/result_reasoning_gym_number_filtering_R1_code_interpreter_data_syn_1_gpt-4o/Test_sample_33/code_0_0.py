# Original list of numbers as strings
numbers = ['-42.2311', '68.942', '-8.309', '10.25', '65.178', '80', '77.083', '-47.192', '25']

# Convert strings to floats and filter numbers greater than -14.5
filtered_numbers = [num for num in numbers if float(num) > -14.5]

# Print the filtered list
print(filtered_numbers)