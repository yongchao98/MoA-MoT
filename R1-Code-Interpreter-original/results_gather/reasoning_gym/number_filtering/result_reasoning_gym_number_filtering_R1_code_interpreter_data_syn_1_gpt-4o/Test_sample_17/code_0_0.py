# Original list of strings
numbers = ['-31', '-35.3010', '2', '35.97']

# Convert strings to floats, filter, and convert back to strings
filtered_numbers = [num for num in numbers if float(num) <= -24]

# Print the result in the required format
print(filtered_numbers)