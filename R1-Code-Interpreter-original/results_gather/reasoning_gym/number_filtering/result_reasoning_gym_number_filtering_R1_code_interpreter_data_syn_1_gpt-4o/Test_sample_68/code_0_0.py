# Original list of strings
original_list = ['-8.1', '41.10', '-13', '7', '47.34']

# Convert strings to floats and filter numbers smaller than -10.1
filtered_list = [str(num) for num in map(float, original_list) if num < -10.1]

# Print the filtered list
print(filtered_list)