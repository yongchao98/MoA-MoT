def find_number_greater_than_threshold(numbers, threshold):
    # Filter numbers greater than the threshold
    filtered_numbers = [num for num in numbers if num > threshold]
    # Return the smallest number from the filtered list
    return min(filtered_numbers) if filtered_numbers else None

# Input list
input_list = [37, 48, 80, 6, 70, 41, 51, 56, 19, 5]
# Threshold
threshold = 50

# Find the number
result = find_number_greater_than_threshold(input_list, threshold)
print([result])