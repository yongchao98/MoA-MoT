# Original list of strings
original_list = ['-12.21', '-3.396', '-23.922', '-33.986', '75.3']

# Convert strings to floats and filter out numbers smaller than 45
filtered_list = [str(num) for num in map(float, original_list) if num >= 45]

# Print the filtered list
print(filtered_list)