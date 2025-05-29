# Original list of numbers as strings
numbers = ['-64.0362', '33.0677', '51', '-84.140', '26.1', '98.11', '88.14']

# Convert strings to floats and filter out numbers larger than -45.7015
filtered_numbers = [num for num in numbers if float(num) <= -45.7015]

# Print the filtered list
print(filtered_numbers)