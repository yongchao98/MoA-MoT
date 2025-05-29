# Original list of numbers as strings
numbers = ['-56.370', '25.1', '22.8962', '-64.54']

# Convert strings to floats and filter out numbers larger than -38.932
filtered_numbers = [num for num in numbers if float(num) <= -38.932]

# Print the filtered list
print(filtered_numbers)