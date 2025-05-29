# Initial positions
positions = {
    "Alice": "Center Midfielder",
    "Bob": "Right Winger",
    "Claire": "Right Midfielder",
    "Dave": "Benchwarmer",
    "Eve": "Cheerleader",
    "Fred": "Goalkeeper",
    "Gertrude": "Left Winger"
}

# List of swaps
swaps = [
    ("Fred", "Bob"),
    ("Bob", "Claire"),
    ("Alice", "Eve"),
    ("Dave", "Claire"),
    ("Fred", "Alice"),
    ("Bob", "Fred"),
    ("Alice", "Gertrude")
]

# Perform swaps
for player1, player2 in swaps:
    positions[player1], positions[player2] = positions[player2], positions[player1]

# Output Gertrude's final position
print(positions["Gertrude"])