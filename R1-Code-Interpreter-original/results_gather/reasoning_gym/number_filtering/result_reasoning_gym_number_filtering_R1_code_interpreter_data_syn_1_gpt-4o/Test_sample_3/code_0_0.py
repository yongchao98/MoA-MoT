# Original list of numbers as strings
numbers = ['50.11', '73.01', '-13.65', '13.3480', '50.98', '30.5795', '86', '-27.1419']

# Convert strings to floats and filter based on the condition
filtered_numbers = [num for num in numbers if float(num) < -25.866]

# Print the filtered list
print(filtered_numbers)