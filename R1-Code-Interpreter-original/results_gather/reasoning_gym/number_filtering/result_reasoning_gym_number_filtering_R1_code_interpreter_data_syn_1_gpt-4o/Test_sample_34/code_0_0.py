# Original list of numbers as strings
numbers = ['-37.6', '57.39', '-48.7824', '-80.962']

# Threshold value
threshold = -55.4733

# Filter and convert the numbers
filtered_numbers = [num for num in numbers if float(num) > threshold]

# Print the filtered list
print(filtered_numbers)