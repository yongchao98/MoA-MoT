input_list = [35, 70, 19, 13, 83, 34, 36, 57, 76]
threshold = 39

# Filter numbers greater than or equal to the threshold
filtered_numbers = [num for num in input_list if num >= threshold]

# Find the smallest number in the filtered list
result = min(filtered_numbers)

print(result)