# Original list of numbers as strings
numbers = ['35.26', '-66.0628', '-63.50', '23.4995', '16.2627']

# Threshold value
threshold = -11.0457

# Filter the list
filtered_numbers = [num for num in numbers if float(num) <= threshold]

# Print the filtered list
print(filtered_numbers)