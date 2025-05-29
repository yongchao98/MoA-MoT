# Original list of numbers as strings
numbers = ['3.2380', '11.6281', '-47', '14.055', '-68', '14.52', '-88', '28.4', '72.5']

# Convert strings to floats, filter, and convert back to strings
filtered_numbers = [num for num in numbers if float(num) >= 9.3]

# Print the filtered list
print(filtered_numbers)