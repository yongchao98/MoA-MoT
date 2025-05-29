# List of numbers to sort
numbers = [26.4, -4.21, 82.0, 95.85, 88.7]

# Sort the numbers in descending order
sorted_numbers = sorted(numbers, reverse=True)

# Convert each number to a string
sorted_numbers_as_strings = [str(num) for num in sorted_numbers]

# Print the result
print(sorted_numbers_as_strings)