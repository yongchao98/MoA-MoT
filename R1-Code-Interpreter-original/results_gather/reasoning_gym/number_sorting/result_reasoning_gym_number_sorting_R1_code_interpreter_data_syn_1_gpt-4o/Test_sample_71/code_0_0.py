# List of numbers
numbers = [86.84, 47.0, 37.82, 47.83, -23.6, 74.0, -2.0, -56.0, -91.0, 38.29]

# Sort the numbers in descending order
sorted_numbers = sorted(numbers, reverse=True)

# Convert each number to a string
sorted_numbers_str = [str(num) for num in sorted_numbers]

# Print the result
print(sorted_numbers_str)