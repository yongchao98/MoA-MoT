# Original list of numbers as strings
numbers = ['-41.184', '-9.99', '99.0', '63.9635']

# Convert strings to floats and filter based on the condition
filtered_numbers = [num for num in numbers if float(num) < -10.831]

# Print the filtered list
print(filtered_numbers)