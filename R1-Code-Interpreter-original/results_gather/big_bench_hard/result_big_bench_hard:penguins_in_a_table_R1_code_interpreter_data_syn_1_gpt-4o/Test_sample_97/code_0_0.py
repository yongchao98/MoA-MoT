# Define the penguins and their weights
penguins = {
    "Louis": 11,
    "Bernard": 13,
    "Vincent": 11,
    "Gwen": 15
}

# Sort the penguins by weight in descending order
sorted_penguins = sorted(penguins.items(), key=lambda x: x[1], reverse=True)

# Get the second heaviest penguin
second_heaviest_penguin = sorted_penguins[1][0]

print(second_heaviest_penguin)