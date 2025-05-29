# Original list of numbers as strings
numbers = ['-45.336', '29', '-28.3276', '-32.32', '31.6609', '77', '11.7', '64', '-17.17']

# Convert strings to floats and filter out numbers smaller than -27.451
filtered_numbers = [num for num in numbers if float(num) >= -27.451]

# Print the filtered list
print(filtered_numbers)