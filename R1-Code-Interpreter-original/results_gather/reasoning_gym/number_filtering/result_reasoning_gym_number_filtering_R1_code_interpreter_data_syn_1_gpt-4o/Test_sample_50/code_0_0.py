# Original list of numbers as strings
numbers = ['53.8791', '-41.6', '58.03', '-69', '-13.3034']

# Convert strings to floats and filter numbers smaller than 37.6
filtered_numbers = [num for num in numbers if float(num) < 37.6]

# Print the filtered list
print(filtered_numbers)