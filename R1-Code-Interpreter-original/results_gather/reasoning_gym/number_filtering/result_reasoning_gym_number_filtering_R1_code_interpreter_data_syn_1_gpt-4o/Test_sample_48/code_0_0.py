# Original list of numbers as strings
numbers = ['-21.212', '-26.1843', '-20', '10.018', '85.7081', '-89.1', '-66.0982']

# Convert strings to floats, filter, and convert back to strings
filtered_numbers = [num for num in numbers if float(num) < 9.3]

# Print the filtered list
print(filtered_numbers)