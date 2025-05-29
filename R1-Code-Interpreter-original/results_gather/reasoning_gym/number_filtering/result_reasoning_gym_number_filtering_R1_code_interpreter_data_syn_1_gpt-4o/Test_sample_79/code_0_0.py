# Original list of numbers as strings
numbers = ['84.23', '-1.216', '-16.6', '25.946', '11.8', '34.086']

# Convert strings to floats, filter, and convert back to strings
filtered_numbers = [str(num) for num in map(float, numbers) if num < 66.23]

# Print the filtered list
print(filtered_numbers)