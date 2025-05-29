# Original list of numbers as strings
numbers = ['-96.6473', '64', '90.66', '60.9', '61.943', '-27.213', '48.07', '38.041', '-44.375', '29.7']

# Convert strings to floats and filter based on the condition
filtered_numbers = [num for num in numbers if float(num) < -45.173]

# Print the filtered list
print(filtered_numbers)