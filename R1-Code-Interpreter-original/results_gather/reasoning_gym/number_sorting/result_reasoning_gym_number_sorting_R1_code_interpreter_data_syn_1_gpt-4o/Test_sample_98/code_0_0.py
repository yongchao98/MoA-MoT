# List of numbers
numbers = [-45.34, 32.0, 29.0, -28.33, -32.3, 31.66, 76.75, -78.9, 90.0]

# Sort the numbers in descending order
sorted_numbers = sorted(numbers, reverse=True)

# Convert each number to a string
sorted_numbers_as_strings = [str(num) for num in sorted_numbers]

# Print the result
print(sorted_numbers_as_strings)