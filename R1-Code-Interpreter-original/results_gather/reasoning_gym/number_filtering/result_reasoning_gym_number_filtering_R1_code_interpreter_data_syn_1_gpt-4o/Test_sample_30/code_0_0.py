# Original list of strings
numbers = ['-87.0', '78.49', '79.55', '-65.8', '-67.69']

# Convert strings to floats and filter numbers greater than 77.1
filtered_numbers = [num for num in numbers if float(num) > 77.1]

# Print the filtered list
print(filtered_numbers)