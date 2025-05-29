# Original list of strings
original_list = ['-46.465', '-47', '-24.1751']

# Convert strings to floats, filter, and convert back to strings
filtered_list = [str(num) for num in map(float, original_list) if num > -39]

# Print the filtered list
print(filtered_list)