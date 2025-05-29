# Original list of numbers as strings
numbers = ['-75.0298', '-19', '-59.30', '40.8', '39.7', '-84.038', '-12', '-83', '85.5', '-79']

# Convert strings to floats and filter numbers greater than 77.03
filtered_numbers = [num for num in numbers if float(num) > 77.03]

# Print the filtered list
print(filtered_numbers)