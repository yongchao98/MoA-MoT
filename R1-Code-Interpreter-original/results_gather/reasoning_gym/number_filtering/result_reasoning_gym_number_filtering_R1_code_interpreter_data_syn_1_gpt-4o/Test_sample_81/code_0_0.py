# Original list of numbers as strings
numbers = ['30.77', '-96.4252', '-30.430']

# Filter numbers smaller than 16
filtered_numbers = [num for num in numbers if float(num) < 16]

# Print the filtered list
print(filtered_numbers)