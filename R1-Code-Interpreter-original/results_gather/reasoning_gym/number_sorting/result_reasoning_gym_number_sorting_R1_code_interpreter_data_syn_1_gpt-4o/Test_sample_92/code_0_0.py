# List of numbers
numbers = [10.0, 89.0, -16.8, -94.0, -20.5, 63.0, 8.61]

# Sort the numbers in descending order
sorted_numbers = sorted(numbers, reverse=True)

# Convert each number to a string
sorted_numbers_str = [str(num) for num in sorted_numbers]

# Print the result
print(sorted_numbers_str)