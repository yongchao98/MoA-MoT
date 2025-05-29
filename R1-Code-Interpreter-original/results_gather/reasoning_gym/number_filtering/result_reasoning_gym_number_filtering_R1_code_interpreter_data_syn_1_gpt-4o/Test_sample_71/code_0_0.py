# Original list of numbers in string format
numbers = ['86.8', '37.8170', '47.829', '-7.77', '74', '-2.2', '-56', '-91', '38.2912', '4.7']

# Convert strings to floats and filter numbers smaller than 41.6152
filtered_numbers = [num for num in numbers if float(num) < 41.6152]

# Print the filtered list
print(filtered_numbers)