# Initial positions
positions = {
    "Alice": "right winger",
    "Bob": "striker",
    "Claire": "benchwarmer",
    "Dave": "fullback",
    "Eve": "left midfielder"
}

# Swaps
swaps = [
    ("Dave", "Bob"),
    ("Claire", "Eve"),
    ("Eve", "Bob"),
    ("Alice", "Claire"),
    ("Alice", "Dave")
]

# Perform swaps
for player1, player2 in swaps:
    positions[player1], positions[player2] = positions[player2], positions[player1]

# Find Claire's final position
claire_final_position = positions["Claire"]
print(claire_final_position)