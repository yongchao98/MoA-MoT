# Original list of numbers as strings
numbers = ['97.6', '-70.01', '68.2', '-39.61']

# Convert strings to floats and filter out numbers smaller than 83.518
filtered_numbers = [num for num in numbers if float(num) >= 83.518]

# Print the filtered list
print(filtered_numbers)