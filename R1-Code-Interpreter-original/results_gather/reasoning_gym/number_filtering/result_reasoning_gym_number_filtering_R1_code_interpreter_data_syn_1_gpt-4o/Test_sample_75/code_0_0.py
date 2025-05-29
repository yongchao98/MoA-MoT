# Original list of numbers as strings
numbers = ['23.962', '-3.50', '91.391', '100', '-49.9769', '-93.05']

# Convert strings to floats and filter numbers greater than 70.1
filtered_numbers = [num for num in numbers if float(num) > 70.1]

# Print the filtered list
print(filtered_numbers)