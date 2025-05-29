# Original list of numbers as strings
numbers = ['19.40', '82', '45.9', '-8.3795', '-67.95', '-10.0']

# Convert strings to floats, filter, and convert back to strings
filtered_numbers = [num for num in numbers if float(num) <= 15.735]

# Print the filtered list
print(filtered_numbers)