# Original list of numbers as strings
numbers = ['-83.078', '-95.45', '42.8746', '-97.885', '-24.3476', '-40.5']

# Convert strings to floats, filter, and convert back to strings
filtered_numbers = [num for num in numbers if float(num) > 4.68]

# Print the filtered list
print(filtered_numbers)