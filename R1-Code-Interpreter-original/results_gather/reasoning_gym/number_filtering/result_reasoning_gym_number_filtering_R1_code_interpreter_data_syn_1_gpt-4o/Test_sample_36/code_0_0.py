# Original list of numbers as strings
numbers = ['88.1608', '-77.8045', '-51.327', '-46.90', '-12.3748', '-27.7493', '30', '-91.643', '-78.79', '-88.9049']

# Convert strings to floats and filter out numbers larger than -59.4
filtered_numbers = [num for num in numbers if float(num) <= -59.4]

# Print the filtered list
print(filtered_numbers)