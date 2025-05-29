# Original list of numbers as strings
numbers = ['-37.067', '85', '64', '40', '-66.72', '-82.61', '-28.6', '1']

# Convert strings to floats and filter out numbers smaller than 58.867
filtered_numbers = [num for num in numbers if float(num) >= 58.867]

# Print the filtered list
print(filtered_numbers)