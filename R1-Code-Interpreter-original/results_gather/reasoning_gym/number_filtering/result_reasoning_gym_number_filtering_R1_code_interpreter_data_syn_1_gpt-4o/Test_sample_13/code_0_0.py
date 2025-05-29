# Original list of numbers as strings
numbers = ['-16.46', '-83.62', '-95.681', '-97.79', '61.53', '71', '80']

# Convert strings to floats and filter out numbers larger than -73.8
filtered_numbers = [num for num in numbers if float(num) <= -73.8]

# Print the filtered list
print(filtered_numbers)