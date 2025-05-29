# Initial positions
positions = {
    "Alice": "right winger",
    "Bob": "left midfielder",
    "Claire": "right midfielder",
    "Dave": "striker",
    "Eve": "fullback"
}

# Swaps
swaps = [
    ("Eve", "Claire"),
    ("Dave", "Bob"),
    ("Eve", "Dave"),
    ("Alice", "Claire"),
    ("Bob", "Eve")
]

# Perform swaps
for player1, player2 in swaps:
    positions[player1], positions[player2] = positions[player2], positions[player1]

# Output Claire's final position
print(positions["Claire"])