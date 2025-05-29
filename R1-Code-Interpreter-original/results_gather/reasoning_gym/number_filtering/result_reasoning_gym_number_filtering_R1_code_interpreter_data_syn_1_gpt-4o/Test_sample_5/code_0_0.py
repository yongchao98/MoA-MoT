# Original list of numbers as strings
numbers = ['87.0627', '-92.91', '-12', '96.502']

# Convert strings to floats and filter out numbers smaller than 20.86
filtered_numbers = [num for num in numbers if float(num) >= 20.86]

# Print the filtered list
print(filtered_numbers)