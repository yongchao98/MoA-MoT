# Initial positions
positions = {
    "Alice": "Right Winger",
    "Bob": "Benchwarmer",
    "Claire": "Right Midfielder",
    "Dave": "Goalkeeper",
    "Eve": "Cheerleader"
}

# Swaps
swaps = [
    ("Bob", "Alice"),
    ("Dave", "Claire"),
    ("Alice", "Claire"),
    ("Eve", "Alice"),
    ("Dave", "Bob")
]

# Perform swaps
for player1, player2 in swaps:
    positions[player1], positions[player2] = positions[player2], positions[player1]

# Output Claire's final position
print(positions["Claire"])