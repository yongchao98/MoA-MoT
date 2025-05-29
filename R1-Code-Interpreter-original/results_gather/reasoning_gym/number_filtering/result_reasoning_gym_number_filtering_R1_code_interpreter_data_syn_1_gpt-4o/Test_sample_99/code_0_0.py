# Original list of numbers as strings
numbers = ['-25.5101', '-82.48', '-97', '-30.3814', '93.1', '-55.5701', '36.583', '95.5']

# Convert strings to floats and filter numbers greater than 85.2175
filtered_numbers = [num for num in numbers if float(num) > 85.2175]

# Print the filtered list
print(filtered_numbers)