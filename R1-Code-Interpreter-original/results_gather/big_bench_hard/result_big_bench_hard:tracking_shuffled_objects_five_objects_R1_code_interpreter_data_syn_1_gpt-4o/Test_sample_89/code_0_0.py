# Initial positions
positions = {
    "Alice": "striker",
    "Bob": "fullback",
    "Claire": "cheerleader",
    "Dave": "benchwarmer",
    "Eve": "goalkeeper"
}

# Swaps
swaps = [
    ("Alice", "Eve"),
    ("Alice", "Claire"),
    ("Dave", "Bob"),
    ("Alice", "Dave"),
    ("Dave", "Eve")
]

# Apply swaps
for player1, player2 in swaps:
    positions[player1], positions[player2] = positions[player2], positions[player1]

# Output Bob's final position
print(positions["Bob"])