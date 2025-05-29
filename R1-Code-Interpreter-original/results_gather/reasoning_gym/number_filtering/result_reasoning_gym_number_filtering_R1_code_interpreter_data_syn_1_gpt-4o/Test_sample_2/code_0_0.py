# Original list of numbers as strings
numbers = ['17.1', '79.729', '-27', '-15.19', '-31.7732', '2.3', '-80.70', '22.1']

# Convert strings to floats and filter out numbers greater than 40
filtered_numbers = [num for num in numbers if float(num) <= 40]

# Print the filtered list
print(filtered_numbers)