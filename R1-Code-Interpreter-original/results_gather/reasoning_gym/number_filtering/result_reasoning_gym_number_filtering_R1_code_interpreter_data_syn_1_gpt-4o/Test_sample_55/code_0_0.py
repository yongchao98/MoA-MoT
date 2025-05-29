# Original list of numbers as strings
numbers = ['47.5830', '-79.56', '28.94', '76', '79']

# Convert strings to floats, filter, and convert back to strings
filtered_numbers = [num for num in numbers if float(num) < -36]

# Print the filtered list
print(filtered_numbers)