# Original list of numbers as strings
numbers = ['26.45', '-4.2066', '81.8', '95.8461', '88.70']

# Convert strings to floats and filter out numbers smaller than 92.4641
filtered_numbers = [num for num in numbers if float(num) >= 92.4641]

# Print the filtered list
print(filtered_numbers)