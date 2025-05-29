# Original list of numbers as strings
numbers = ['58.9297', '-50.447', '42.7027', '37.7693', '-16.65', '49', '63', '-70', '5.731']

# Convert strings to floats and filter out numbers smaller than 56.3
filtered_numbers = [num for num in numbers if float(num) >= 56.3]

# Print the filtered list
print(filtered_numbers)