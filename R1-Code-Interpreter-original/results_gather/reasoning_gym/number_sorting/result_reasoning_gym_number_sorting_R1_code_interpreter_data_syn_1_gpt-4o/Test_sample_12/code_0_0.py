# List of numbers
numbers = [4.0, 72.0, -24.1, -94.0, 14.0, -68.66, 37.8, 38.7, 18.2]

# Sort the numbers in descending order
sorted_numbers = sorted(numbers, reverse=True)

# Convert each number to a string
sorted_numbers_str = [str(num) for num in sorted_numbers]

# Print the result
print(sorted_numbers_str)