# Original list of numbers as strings
numbers = ['56.4', '-31.2032', '23.96', '47.9205']

# Convert strings to floats and filter out numbers larger than -19.90
filtered_numbers = [num for num in numbers if float(num) <= -19.90]

# Print the filtered list
print(filtered_numbers)