# Initial positions
positions = {
    "Alice": "center midfielder",
    "Bob": "fullback",
    "Claire": "goalkeeper",
    "Dave": "left winger",
    "Eve": "right winger"
}

# Swaps
swaps = [
    ("Eve", "Claire"),
    ("Alice", "Dave"),
    ("Claire", "Bob"),
    ("Eve", "Alice"),
    ("Claire", "Eve")
]

# Perform swaps
for player1, player2 in swaps:
    positions[player1], positions[player2] = positions[player2], positions[player1]

# Output Claire's final position
print(positions["Claire"])