# List of numbers as strings
numbers = [
    "585,747,269.003",
    "585747275.870000",
    "585,747,276.9760001",
    "585,747,276.6470001"
]

# Convert the numbers to float after removing commas
float_numbers = [float(num.replace(',', '')) for num in numbers]

# Find the smallest number
smallest_number = min(float_numbers)

# Print the smallest number
print(smallest_number)