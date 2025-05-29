# Original list of numbers as strings
numbers = ['-11.15', '-3', '-23', '-84.9', '-41.85', '68.5067', '-88.0', '-34.37', '20', '94.5472']

# Convert strings to floats and filter out numbers smaller than 42.7
filtered_numbers = [num for num in numbers if float(num) >= 42.7]

# Print the filtered list
print(filtered_numbers)