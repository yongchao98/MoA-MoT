# Original list of strings
original_list = ['74.4', '-22', '43']

# Convert strings to floats and filter numbers greater than 66
filtered_list = [str(num) for num in map(float, original_list) if num > 66]

# Print the filtered list
print(filtered_list)