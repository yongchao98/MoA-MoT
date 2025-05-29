# Original list of numbers as strings
numbers = ['42.719', '31.8960', '-38.04', '60.9952', '77.67', '-47.05']

# Convert strings to floats and filter numbers greater than -0.9340
filtered_numbers = [num for num in numbers if float(num) > -0.9340]

# Print the filtered list
print(filtered_numbers)