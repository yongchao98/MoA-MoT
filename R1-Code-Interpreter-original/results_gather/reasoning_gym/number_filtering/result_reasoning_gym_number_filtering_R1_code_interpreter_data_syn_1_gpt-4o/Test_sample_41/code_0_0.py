# Original list of numbers as strings
numbers = ['-75.558', '12.830', '-40', '-2.5532', '0.101', '-77.92', '-18']

# Threshold for filtering
threshold = -40.2034

# Convert strings to floats, filter, and convert back to strings
filtered_numbers = [num for num in numbers if float(num) >= threshold]

# Print the result
print(filtered_numbers)