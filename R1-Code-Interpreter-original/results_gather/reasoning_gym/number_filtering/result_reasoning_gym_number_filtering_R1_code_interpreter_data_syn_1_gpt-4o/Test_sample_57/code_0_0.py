# Original list of numbers as strings
numbers = ['53.57', '-72', '-70.085', '-33.1']

# Convert strings to floats and filter numbers greater than -34.37
filtered_numbers = [num for num in numbers if float(num) > -34.37]

# Print the filtered list
print(filtered_numbers)