# Original list of numbers as strings
numbers = ['-76.781', '94.419', '-46.6300', '-14.1067']

# Convert strings to floats and filter numbers larger than -29.2
filtered_numbers = [num for num in numbers if float(num) > -29.2]

# Print the filtered list
print(filtered_numbers)