# Original list of strings
original_list = ['4', '-64.7', '-42.1', '-77', '-79.9640', '37.76', '38.702', '18.20', '-28.34']

# Convert strings to floats and filter numbers greater than 19.8962
filtered_list = [str(num) for num in map(float, original_list) if num > 19.8962]

# Print the filtered list
print(filtered_list)