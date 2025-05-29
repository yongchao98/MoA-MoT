# Original list of numbers as strings
numbers = ['-36.15', '-63.6', '-29.718', '60.3380']

# Convert strings to floats and filter out numbers greater than 50.936
filtered_numbers = [num for num in numbers if float(num) <= 50.936]

# Print the filtered list
print(filtered_numbers)