# Define the penguins with their respective heights
penguins = {
    "Louis": 50,
    "Bernard": 80,
    "Vincent": 60,
    "Gwen": 70
}

# Find the penguin with the minimum height
shortest_penguin = min(penguins, key=penguins.get)

# Output the shortest penguin
print(shortest_penguin)