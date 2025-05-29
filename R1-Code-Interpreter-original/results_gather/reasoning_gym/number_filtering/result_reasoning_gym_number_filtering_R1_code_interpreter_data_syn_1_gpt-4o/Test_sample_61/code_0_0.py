# Original list of numbers as strings
numbers = ['84.00', '77', '-28.9', '-21.7', '16.3323', '11', '35.75', '-20', '31.1', '50.391']

# Convert strings to floats, filter, and convert back to strings
filtered_numbers = [str(num) for num in map(float, numbers) if num >= -10.61]

# Print the filtered list
print(filtered_numbers)