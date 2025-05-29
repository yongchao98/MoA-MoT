# Define the penguins and their ages
penguins = {
    "Louis": 7,
    "Bernard": 5,
    "Vincent": 9,
    "Gwen": 8
}

# Sort the penguins by age
sorted_penguins = sorted(penguins.items(), key=lambda x: x[1])

# Get the second youngest penguin
second_youngest = sorted_penguins[1][0]

# Print the result
print(second_youngest)