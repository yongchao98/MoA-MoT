# Original list of numbers as strings
numbers = ['40', '-35.8', '-55.50', '16.444', '7.2762', '-55.0834']

# Convert strings to floats, filter, and convert back to strings
filtered_numbers = [num for num in numbers if float(num) < 1.924]

# Print the filtered list
print(filtered_numbers)