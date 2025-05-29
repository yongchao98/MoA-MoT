# Original list of strings
original_list = ['-95.00', '-51.0', '47.2942', '-82.612']

# Convert strings to floats and filter numbers greater than -90
filtered_list = [num for num in original_list if float(num) > -90]

# Print the filtered list
print(filtered_list)