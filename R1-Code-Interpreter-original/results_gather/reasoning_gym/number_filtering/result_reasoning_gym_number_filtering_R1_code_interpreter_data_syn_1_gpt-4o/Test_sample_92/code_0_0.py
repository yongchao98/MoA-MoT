# Original list of numbers as strings
numbers = ['9.7', '89', '-16.76', '-93.7', '-20.51', '63', '8.6111']

# Convert strings to floats and filter numbers larger than -15.9
filtered_numbers = [num for num in numbers if float(num) > -15.9]

# Print the filtered list
print(filtered_numbers)