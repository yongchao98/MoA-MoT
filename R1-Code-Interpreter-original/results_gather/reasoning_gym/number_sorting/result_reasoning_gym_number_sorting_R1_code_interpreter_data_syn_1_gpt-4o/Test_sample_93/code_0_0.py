# List of numbers to sort
numbers = [-55.49, -39.81, -64.0, 21.5, 26.0, 42.0]

# Sort the numbers in descending order
sorted_numbers = sorted(numbers, reverse=True)

# Convert each number to a string
sorted_numbers_as_strings = [str(num) for num in sorted_numbers]

# Print the result
print(sorted_numbers_as_strings)