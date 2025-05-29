# List of numbers to be sorted
numbers = [-58.86, 96.0, -91.8, -18.9, 61.1, 93.35]

# Sorting the numbers in descending order
sorted_numbers = sorted(numbers, reverse=True)

# Converting each number to a string
sorted_numbers_as_strings = [str(num) for num in sorted_numbers]

# Output the result
print(sorted_numbers_as_strings)