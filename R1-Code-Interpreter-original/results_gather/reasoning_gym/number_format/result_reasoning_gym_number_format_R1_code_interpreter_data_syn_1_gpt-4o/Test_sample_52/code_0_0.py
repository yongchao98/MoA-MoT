# List of numbers as strings to handle commas
numbers = [
    "914,061,876.456",
    "914061867.203000",
    "914,061,866.544",
    "914,061,886.3249999"
]

# Remove commas and convert to float
numbers_float = [float(num.replace(',', '')) for num in numbers]

# Find the smallest number
smallest_number = min(numbers_float)

# Print the smallest number
print(smallest_number)