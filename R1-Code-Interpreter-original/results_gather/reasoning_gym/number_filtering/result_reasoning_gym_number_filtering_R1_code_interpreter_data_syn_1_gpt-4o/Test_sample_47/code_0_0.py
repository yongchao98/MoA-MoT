# Original list of numbers as strings
numbers = ['-12.4993', '-67.4', '44.744', '-69.00', '-94']

# Convert strings to floats and filter numbers smaller than 41.193
filtered_numbers = [num for num in numbers if float(num) < 41.193]

# Print the filtered list
print(filtered_numbers)