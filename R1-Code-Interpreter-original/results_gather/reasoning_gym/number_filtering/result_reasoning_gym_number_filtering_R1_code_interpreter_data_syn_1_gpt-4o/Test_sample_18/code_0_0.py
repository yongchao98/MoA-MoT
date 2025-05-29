# Original list of numbers as strings
numbers = ['71.14', '27.5', '38.51', '-82.96', '-55.1214', '-83.0', '-31', '95.94', '-56', '94.48']

# Threshold value
threshold = 93.8156

# Filter numbers smaller than the threshold
filtered_numbers = [num for num in numbers if float(num) < threshold]

# Print the filtered list
print(filtered_numbers)