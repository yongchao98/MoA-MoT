# Original list of numbers as strings
numbers = ['-32.3824', '-23', '-68.118', '34.5673', '-99.71', '-42.26', '-55.7', '20.7', '2.79', '75.0443']

# Convert strings to floats and filter out numbers greater than 68.24
filtered_numbers = [num for num in numbers if float(num) <= 68.24]

# Print the filtered list
print(filtered_numbers)