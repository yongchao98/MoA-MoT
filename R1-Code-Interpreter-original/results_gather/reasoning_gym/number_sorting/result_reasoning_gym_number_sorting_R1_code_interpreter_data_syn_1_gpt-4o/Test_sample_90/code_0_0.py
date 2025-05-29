# List of numbers to sort
numbers = [19.4, 82.0, 46.0, -8.38, 46.77, -67.9]

# Sort the numbers in descending order
sorted_numbers = sorted(numbers, reverse=True)

# Convert each number to a string
sorted_numbers_as_strings = [str(num) for num in sorted_numbers]

# Print the result
print(sorted_numbers_as_strings)