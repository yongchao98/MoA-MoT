# List of numbers
numbers = [17.0, -91.67, 34.69, 85.22, 35.87, 64.0, 79.5, 85.4, -45.08, 72.3]

# Sort the numbers in descending order
sorted_numbers = sorted(numbers, reverse=True)

# Convert each number to a string
sorted_numbers_str = [str(num) for num in sorted_numbers]

# Print the result
print(sorted_numbers_str)