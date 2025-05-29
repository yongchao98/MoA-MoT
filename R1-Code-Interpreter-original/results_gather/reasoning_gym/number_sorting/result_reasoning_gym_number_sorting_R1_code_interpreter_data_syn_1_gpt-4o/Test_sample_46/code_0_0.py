# List of numbers to sort
numbers = [-80.8, 46.1, 49.0, 31.1, -19.86, -62.26]

# Sort the numbers in descending order
sorted_numbers = sorted(numbers, reverse=True)

# Convert each number to a string
sorted_numbers_str = [str(num) for num in sorted_numbers]

# Print the result
print(sorted_numbers_str)