# List of numbers to sort
numbers = [-96.65, 64.0, 90.7, 61.0, 61.9, -27.2, 48.1, 38.04, -6.37, -3.0]

# Sort the numbers in ascending order
sorted_numbers = sorted(numbers)

# Convert each number to a string
sorted_strings = [str(num) for num in sorted_numbers]

# Print the result in the specified format
print(f"<<<{sorted_strings}>>>")