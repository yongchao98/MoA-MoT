# Initial positions
positions = {
    "Alice": "goalkeeper",
    "Bob": "left winger",
    "Claire": "benchwarmer",
    "Dave": "cheerleader",
    "Eve": "left midfielder"
}

# Swaps
swaps = [
    ("Eve", "Alice"),
    ("Bob", "Eve"),
    ("Alice", "Dave"),
    ("Bob", "Eve"),
    ("Eve", "Claire")
]

# Apply swaps
for player1, player2 in swaps:
    positions[player1], positions[player2] = positions[player2], positions[player1]

# Output Bob's final position
print(positions["Bob"])