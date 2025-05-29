# Original list of numbers as strings
numbers = ['-54.197', '9', '-54.468', '91.2740', '83', '57.986', '1.197']

# Convert strings to floats and filter out numbers smaller than 67.8
filtered_numbers = [num for num in numbers if float(num) >= 67.8]

# Print the filtered list
print(filtered_numbers)