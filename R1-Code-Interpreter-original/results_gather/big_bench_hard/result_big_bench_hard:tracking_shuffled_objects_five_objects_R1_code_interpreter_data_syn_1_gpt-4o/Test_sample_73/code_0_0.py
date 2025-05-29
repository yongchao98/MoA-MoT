# Initial positions
positions = {
    "Alice": "Striker",
    "Bob": "Right Winger",
    "Claire": "Goalkeeper",
    "Dave": "Left Winger",
    "Eve": "Left Midfielder"
}

# Swaps
swaps = [
    ("Dave", "Bob"),
    ("Alice", "Bob"),
    ("Eve", "Dave"),
    ("Bob", "Dave"),
    ("Claire", "Alice")
]

# Apply swaps
for player1, player2 in swaps:
    positions[player1], positions[player2] = positions[player2], positions[player1]

# Final position of Claire
print(positions["Claire"])