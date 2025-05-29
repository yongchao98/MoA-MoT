# Define the penguin data as a list of strings, each representing a line in the table
penguin_data = [
    "name, age, height (cm), weight (kg)",
    "Louis, 7, 50, 11",
    "Bernard, 5, 80, 13",
    "Vincent, 9, 60, 11",
    "Gwen, 8, 70, 15"
]

# Count the number of penguins by subtracting 1 for the header
number_of_penguins = len(penguin_data) - 1

# Output the result
print(number_of_penguins)