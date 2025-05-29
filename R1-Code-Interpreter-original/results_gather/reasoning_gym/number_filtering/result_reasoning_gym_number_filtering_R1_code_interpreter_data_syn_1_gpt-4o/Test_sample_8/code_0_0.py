# Original list of numbers as strings
numbers = ['15.8975', '-93.63', '32.9', '-74.33', '-44.576', '63.7749', '79.2', '-88.1', '23.916', '70']

# Convert strings to floats and filter
filtered_numbers = [num for num in numbers if float(num) < -63.1164]

# Print the filtered list
print(filtered_numbers)