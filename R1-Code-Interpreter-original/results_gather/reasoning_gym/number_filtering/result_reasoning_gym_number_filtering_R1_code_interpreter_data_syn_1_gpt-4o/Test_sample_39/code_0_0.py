# Original list of strings
original_list = ['2', '-47.0', '64.4', '-47.3', '2.4465', '-61', '62.5076', '-26.866']

# Convert strings to floats and filter numbers greater than 34.39
filtered_list = [str(num) for num in map(float, original_list) if num > 34.39]

# Print the filtered list
print(filtered_list)