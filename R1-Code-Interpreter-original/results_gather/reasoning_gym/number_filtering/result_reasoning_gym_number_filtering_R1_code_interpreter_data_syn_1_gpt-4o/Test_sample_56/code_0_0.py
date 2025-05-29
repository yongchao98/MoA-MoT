# Original list of numbers as strings
numbers = ['94.74', '80', '5', '70.496', '-73', '37.6622', '-65.336', '-35.264', '-72.1']

# Convert strings to floats and filter numbers larger than -52.89
filtered_numbers = [num for num in numbers if float(num) > -52.89]

# Print the filtered list
print(filtered_numbers)