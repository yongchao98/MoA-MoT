# Initial positions
positions = {
    "Alice": "cheerleader",
    "Bob": "left winger",
    "Claire": "goalkeeper",
    "Dave": "right midfielder",
    "Eve": "center midfielder",
    "Fred": "benchwarmer",
    "Gertrude": "striker"
}

# Swaps
swaps = [
    ("Fred", "Claire"),
    ("Gertrude", "Alice"),
    ("Fred", "Dave"),
    ("Fred", "Claire"),
    ("Alice", "Bob"),
    ("Dave", "Bob"),
    ("Fred", "Eve")
]

# Apply swaps
for player1, player2 in swaps:
    positions[player1], positions[player2] = positions[player2], positions[player1]

# Output Bob's final position
print(positions["Bob"])