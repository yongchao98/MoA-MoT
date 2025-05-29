# Original list of numbers as strings
numbers = ['-55.4935', '-39.8', '-29.7863', '-39', '-84.9', '-6.093']

# Convert strings to floats, filter, and convert back to strings
filtered_numbers = [num for num in numbers if float(num) < -34.80]

# Print the filtered list
print(filtered_numbers)