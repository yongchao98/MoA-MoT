# Original list of strings
original_list = ['-20.0011', '17.4', '26.9', '17']

# Convert strings to floats and filter
filtered_list = [str(num) for num in map(float, original_list) if num > -18.58]

# Print the filtered list
print(filtered_list)