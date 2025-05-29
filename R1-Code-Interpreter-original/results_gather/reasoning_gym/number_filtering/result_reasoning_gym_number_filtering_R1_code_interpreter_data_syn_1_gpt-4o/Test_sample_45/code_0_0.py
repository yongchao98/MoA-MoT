# Original list of numbers as strings
numbers = ['-34.8', '-61', '-41.5350', '99.9', '-70.9057', '-99.51', '-94.330']

# Convert strings to floats, filter, and convert back to strings
filtered_numbers = [num for num in numbers if float(num) < -65.1]

# Print the filtered list
print(filtered_numbers)