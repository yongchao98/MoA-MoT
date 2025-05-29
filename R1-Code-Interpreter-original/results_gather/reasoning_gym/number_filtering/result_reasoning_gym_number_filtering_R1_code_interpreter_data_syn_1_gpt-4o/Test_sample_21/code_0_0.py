# Original list of numbers as strings
numbers = ['-8.703', '87.3963', '-27.859', '46', '82.26', '49.36']

# Convert strings to floats and filter numbers greater than 74
filtered_numbers = [num for num in numbers if float(num) > 74]

# Print the filtered list
print(filtered_numbers)