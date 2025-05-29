# Initial positions
positions = {
    "Alice": "fullback",
    "Bob": "right winger",
    "Claire": "center midfielder",
    "Dave": "striker",
    "Eve": "benchwarmer"
}

# Swaps
swaps = [
    ("Alice", "Bob"),
    ("Dave", "Bob"),
    ("Alice", "Claire"),
    ("Bob", "Claire"),
    ("Alice", "Eve")
]

# Apply swaps
for player1, player2 in swaps:
    positions[player1], positions[player2] = positions[player2], positions[player1]

# Output Bob's final position
print(positions["Bob"])