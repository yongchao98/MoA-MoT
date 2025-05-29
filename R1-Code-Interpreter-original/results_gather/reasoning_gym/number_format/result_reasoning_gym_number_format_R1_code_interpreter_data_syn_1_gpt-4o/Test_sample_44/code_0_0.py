# List of numbers as strings to handle commas
numbers = [
    "463258311.082000",
    "463,258,308.875",
    "463,258,319.254",
    "463,258,302.137"
]

# Remove commas and convert to float
numbers = [float(num.replace(',', '')) for num in numbers]

# Find the smallest number
smallest_number = min(numbers)

# Print the smallest number
print(smallest_number)