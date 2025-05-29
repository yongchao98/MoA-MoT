# Original list of numbers as strings
numbers = ['16', '90.32', '4.8074', '-33.1', '70.0', '50.2', '-80.002', '-29.7']

# Convert strings to floats and filter out numbers smaller than 64.2251
filtered_numbers = [num for num in numbers if float(num) >= 64.2251]

# Print the filtered list
print(filtered_numbers)