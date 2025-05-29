# Original list of numbers as strings
numbers = ['-88', '91.85', '75', '0.5', '-50.471', '-27.10', '8.5550', '-18.430']

# Convert strings to floats, filter, and convert back to strings
filtered_numbers = [num for num in numbers if float(num) < 88.53]

# Print the filtered list
print(filtered_numbers)