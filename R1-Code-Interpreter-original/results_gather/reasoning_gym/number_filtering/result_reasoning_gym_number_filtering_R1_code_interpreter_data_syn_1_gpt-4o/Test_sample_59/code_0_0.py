# Original list of numbers as strings
numbers = ['18.1', '-19.550', '68.866', '-6.91']

# Convert strings to floats and filter out numbers smaller than 1.9946
filtered_numbers = [num for num in numbers if float(num) >= 1.9946]

# Print the filtered list
print(filtered_numbers)