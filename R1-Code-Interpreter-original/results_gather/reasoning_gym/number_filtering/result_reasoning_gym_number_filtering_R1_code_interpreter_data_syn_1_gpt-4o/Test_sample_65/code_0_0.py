# Original list of numbers as strings
numbers = ['-14.21', '83', '14.8696', '55', '77.65', '-26']

# Convert strings to floats and filter out numbers greater than 43.398
filtered_numbers = [num for num in numbers if float(num) <= 43.398]

# Print the filtered list
print(filtered_numbers)