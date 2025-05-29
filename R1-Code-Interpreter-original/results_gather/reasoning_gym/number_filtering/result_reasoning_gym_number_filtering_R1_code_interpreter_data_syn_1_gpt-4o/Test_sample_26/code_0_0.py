# Original list of numbers as strings
numbers = ['-58.9', '-91.759', '-18.92', '61.08', '93.3462', '-16.1939']

# Convert strings to floats and filter out numbers larger than -72.395
filtered_numbers = [num for num in numbers if float(num) <= -72.395]

# Print the filtered list
print(filtered_numbers)