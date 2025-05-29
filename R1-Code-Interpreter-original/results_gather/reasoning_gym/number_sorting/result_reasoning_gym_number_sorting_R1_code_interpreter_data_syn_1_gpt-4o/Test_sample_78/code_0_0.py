# List of numbers
numbers = [58.93, -50.4, 42.7, 8.0, 8.86, -37.0, -76.74, 82.0, -12.0]

# Sort the numbers in descending order
sorted_numbers = sorted(numbers, reverse=True)

# Convert each number to a string
sorted_strings = [str(num) for num in sorted_numbers]

# Print the result
print(sorted_strings)