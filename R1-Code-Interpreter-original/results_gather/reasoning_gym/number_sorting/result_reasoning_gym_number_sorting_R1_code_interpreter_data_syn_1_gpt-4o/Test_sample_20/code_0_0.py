# List of numbers to sort
numbers = [-89.1, 2.5, -19.0, -68.2, 66.1, -12.0, -93.2]

# Sort the numbers in descending order
sorted_numbers = sorted(numbers, reverse=True)

# Convert each number to a string
sorted_numbers_as_strings = [str(num) for num in sorted_numbers]

# Print the result
print(sorted_numbers_as_strings)