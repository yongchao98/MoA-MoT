# Original list of strings
numbers = ['-61', '7.47', '50.8608', '98.29', '37.463', '-87.6494', '-71']

# Convert strings to floats, filter, and convert back to strings
filtered_numbers = [num for num in numbers if float(num) < -34.7]

# Print the result in the required format
print(f"<<<{filtered_numbers}>>>")