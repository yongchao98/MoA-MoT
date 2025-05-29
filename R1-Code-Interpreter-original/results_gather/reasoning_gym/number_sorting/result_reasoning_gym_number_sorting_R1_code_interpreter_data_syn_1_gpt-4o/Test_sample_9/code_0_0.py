# List of numbers to be sorted
numbers = [-34.0, 80.4, 15.59, -43.3, -23.29, 53.0, 56.24, 100.0, -12.0]

# Step 1: Sort the numbers in ascending order
sorted_numbers = sorted(numbers)

# Step 2: Convert each number to a string
sorted_strings = [str(num) for num in sorted_numbers]

# Step 3: Print the formatted output
print(sorted_strings)