# Original list of numbers as strings
numbers = ['15.613', '-50.5', '-5.631']

# Convert strings to floats and filter based on the condition
filtered_numbers = [num for num in numbers if float(num) < -35.4]

# Print the filtered list
print(filtered_numbers)