# Original list of numbers as strings
numbers = ['-63.7', '-20.0', '12', '43.1867', '-33.793', '62']

# Convert strings to floats and filter based on the condition
filtered_numbers = [num for num in numbers if float(num) < -47.60]

# Print the filtered list
print(filtered_numbers)