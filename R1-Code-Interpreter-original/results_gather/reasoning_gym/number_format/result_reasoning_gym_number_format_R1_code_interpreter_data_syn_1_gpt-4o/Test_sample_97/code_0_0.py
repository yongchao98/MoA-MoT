# List of numbers as strings
numbers = [
    "46174406.296000",
    "46,174,414.544999994",
    "46174410.191000",
    "46174404.272000"
]

# Remove commas and convert to float
numbers = [float(num.replace(',', '')) for num in numbers]

# Find the largest number
largest_number = max(numbers)

# Print the largest number
print(largest_number)