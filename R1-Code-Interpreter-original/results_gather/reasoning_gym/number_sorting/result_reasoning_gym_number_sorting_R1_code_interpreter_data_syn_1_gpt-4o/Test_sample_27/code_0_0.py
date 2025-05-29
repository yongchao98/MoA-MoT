# List of numbers to sort
numbers = [-83.08, 89.0, -72.8, 42.87, 24.3, -24.35]

# Sort the numbers in descending order
sorted_numbers = sorted(numbers, reverse=True)

# Convert each number to a string
sorted_numbers_as_strings = [str(num) for num in sorted_numbers]

# Print the result
print(sorted_numbers_as_strings)