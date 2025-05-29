# Original list of strings
original_list = ['-43.3', '-46.876', '-6.56', '95', '58', '33.2', '-69']

# Convert strings to floats and filter numbers smaller than -30
filtered_list = [str(num) for num in map(float, original_list) if num < -30]

# Print the filtered list
print(filtered_list)