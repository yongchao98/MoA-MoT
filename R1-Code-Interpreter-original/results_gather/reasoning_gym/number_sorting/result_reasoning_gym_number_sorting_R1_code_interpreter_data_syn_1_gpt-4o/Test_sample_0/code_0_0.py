# List of numbers to sort
numbers = [85.0, -39.0, -0.75, -92.3]

# Sort the numbers in descending order
sorted_numbers = sorted(numbers, reverse=True)

# Convert each number to a string
sorted_numbers_str = [str(num) for num in sorted_numbers]

# Print the result
print(sorted_numbers_str)