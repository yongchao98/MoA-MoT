# Initial positions
positions = {
    "Alice": "left midfielder",
    "Bob": "cheerleader",
    "Claire": "right midfielder",
    "Dave": "center midfielder",
    "Eve": "goalkeeper",
    "Fred": "striker",
    "Gertrude": "benchwarmer"
}

# List of swaps
swaps = [
    ("Dave", "Bob"),
    ("Alice", "Eve"),
    ("Claire", "Gertrude"),
    ("Fred", "Alice"),
    ("Bob", "Claire"),
    ("Dave", "Eve"),
    ("Bob", "Alice")
]

# Perform swaps
for player1, player2 in swaps:
    positions[player1], positions[player2] = positions[player2], positions[player1]

# Output Fred's final position
print(positions["Fred"])