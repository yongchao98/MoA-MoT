# Original list of numbers as strings
numbers = ['-14', '-86.754', '-25.0137', '33.0071', '-31.09', '42.0210', '45.319', '98.0809', '-19.8521']

# Convert strings to floats and filter out numbers larger than -71.160
filtered_numbers = [num for num in numbers if float(num) <= -71.160]

# Print the filtered list
print(filtered_numbers)