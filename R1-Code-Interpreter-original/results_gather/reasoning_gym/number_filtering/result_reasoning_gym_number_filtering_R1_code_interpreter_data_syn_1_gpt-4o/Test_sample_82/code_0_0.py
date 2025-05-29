# Original list of numbers as strings
numbers = ['12.54', '33.2', '25', '-69.2921', '34.95', '16.3289']

# Convert strings to floats and filter numbers less than 17.98
filtered_numbers = [num for num in numbers if float(num) < 17.98]

# Print the filtered list
print(filtered_numbers)