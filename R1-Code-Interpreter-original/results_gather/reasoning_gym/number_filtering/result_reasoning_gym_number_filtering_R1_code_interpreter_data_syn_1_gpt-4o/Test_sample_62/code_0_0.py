# Original list of numbers as strings
numbers = ['-75.18', '-40', '37.46', '90.406', '19.02']

# Convert strings to floats, filter, and convert back to strings
filtered_numbers = [num for num in numbers if float(num) < 84]

# Print the result
print(filtered_numbers)