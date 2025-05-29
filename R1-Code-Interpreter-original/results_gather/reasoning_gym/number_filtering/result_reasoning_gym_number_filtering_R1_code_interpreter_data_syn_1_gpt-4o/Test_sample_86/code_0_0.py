# Original list of numbers as strings
numbers = ['-37.57', '-41', '-39', '21.5622', '93.58']

# Convert strings to floats and filter numbers greater than 35.6166
filtered_numbers = [num for num in numbers if float(num) > 35.6166]

# Print the filtered list
print(filtered_numbers)