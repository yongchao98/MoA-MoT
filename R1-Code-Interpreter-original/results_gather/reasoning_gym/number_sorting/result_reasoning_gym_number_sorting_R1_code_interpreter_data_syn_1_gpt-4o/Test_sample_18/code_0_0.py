# List of numbers
numbers = [71.1, 27.0, 38.5, -83.0, -55.12, -83.0, -31.0, 95.9, -56.0, 94.5]

# Sort the numbers in descending order
sorted_numbers = sorted(numbers, reverse=True)

# Convert each number to a string
sorted_numbers_as_strings = [str(num) for num in sorted_numbers]

# Print the result
print(sorted_numbers_as_strings)