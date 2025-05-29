# Original list of string numbers
numbers = ['-26.4496', '-96.1378', '-35.756']

# Convert strings to floats and filter
filtered_numbers = [num for num in numbers if float(num) < -75.911]

# Print the filtered list
print(filtered_numbers)