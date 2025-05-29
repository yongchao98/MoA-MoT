# Original list of numbers as strings
numbers = ['-33.6', '80.401', '15.59', '10.612', '70.3927', '-98.1', '32', '-12.1', '97.33']

# Convert strings to floats and filter numbers larger than -65.3634
filtered_numbers = [num for num in numbers if float(num) > -65.3634]

# Print the filtered list
print(filtered_numbers)