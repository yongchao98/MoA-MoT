# Define the table as a list of strings, each representing a line
table = [
    "name, age, height (cm), weight (kg)",
    "Louis, 7, 50, 11",
    "Bernard, 5, 80, 13",
    "Vincent, 9, 60, 11",
    "Gwen, 8, 70, 15"
]

# The number of penguins is the total number of lines minus the header
number_of_penguins = len(table) - 1

# Output the result
print(number_of_penguins)