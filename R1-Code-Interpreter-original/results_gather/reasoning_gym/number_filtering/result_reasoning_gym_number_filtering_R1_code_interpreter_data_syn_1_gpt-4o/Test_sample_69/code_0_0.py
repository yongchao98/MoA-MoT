# Original list of numbers as strings
numbers = ['81.6756', '84.795', '-90.2958', '-57.1', '83.507', '38.5']

# Convert strings to floats and filter out numbers larger than -32.684
filtered_numbers = [num for num in numbers if float(num) <= -32.684]

# Print the filtered list
print(filtered_numbers)