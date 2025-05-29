# Original list of strings
original_list = ['85.3', '-39.4', '-1', '-79.9529']

# Convert strings to floats and filter out numbers smaller than -25
filtered_list = [str(num) for num in map(float, original_list) if num >= -25]

# Print the filtered list
print(filtered_list)