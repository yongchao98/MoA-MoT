# Original list of numbers as strings
numbers = ['-80.8', '-86.66', '-34.3928']

# Convert strings to floats and filter out numbers smaller than -64.994
filtered_numbers = [num for num in numbers if float(num) >= -64.994]

# Print the filtered list
print(filtered_numbers)