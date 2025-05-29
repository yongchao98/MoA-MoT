# Original list of numbers as strings
numbers = ['-73.6217', '-39.5500', '-61.456', '-67', '-1.4', '-56.3901', '11.3127', '-69.0']

# Convert strings to floats and filter numbers larger than -71.248
filtered_numbers = [num for num in numbers if float(num) > -71.248]

# Print the filtered list
print(filtered_numbers)