# Original list of numbers as strings
numbers = ['-7.348', '-59.99', '-89.4', '-43.674', '-2', '-33.07', '-54.5153', '-87.918']

# Threshold value
threshold = -37.8206

# Filter the list
filtered_numbers = [num for num in numbers if float(num) < threshold]

# Print the filtered list
print(filtered_numbers)