# Original list of numbers as strings
numbers = ['67.6981', '-67.6', '66.52', '85.271', '9.658', '46.2944']

# Convert strings to floats, filter, and convert back to strings
filtered_numbers = [num for num in numbers if float(num) >= 12.43]

# Print the filtered list
print(filtered_numbers)