# Original list of numbers as strings
numbers = ['-87.4950', '-9.27', '-48.514', '-22', '59', '91', '-47.239', '45']

# Convert strings to floats and filter numbers smaller than 76.26
filtered_numbers = [num for num in numbers if float(num) < 76.26]

# Print the filtered list
print(filtered_numbers)