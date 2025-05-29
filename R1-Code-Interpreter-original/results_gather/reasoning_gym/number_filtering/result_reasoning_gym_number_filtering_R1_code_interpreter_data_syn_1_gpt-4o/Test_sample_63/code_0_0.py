# Original list of strings
original_list = ['89', '-72', '20.23', '61', '13.6020']

# Convert strings to floats, filter, and convert back to strings
filtered_list = [str(num) for num in map(float, original_list) if num >= 26]

# Print the filtered list
print(filtered_list)