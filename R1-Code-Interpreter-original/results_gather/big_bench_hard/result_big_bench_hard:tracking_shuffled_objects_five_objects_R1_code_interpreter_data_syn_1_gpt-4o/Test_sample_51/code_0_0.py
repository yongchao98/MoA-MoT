# Initial positions
positions = {
    "Alice": "Left Winger",
    "Bob": "Fullback",
    "Claire": "Benchwarmer",
    "Dave": "Goalkeeper",
    "Eve": "Left Midfielder"
}

# Swaps
swaps = [
    ("Bob", "Claire"),
    ("Bob", "Dave"),
    ("Eve", "Dave"),
    ("Dave", "Claire"),
    ("Alice", "Dave")
]

# Apply swaps
for player1, player2 in swaps:
    positions[player1], positions[player2] = positions[player2], positions[player1]

# Final position of Claire
print(positions["Claire"])