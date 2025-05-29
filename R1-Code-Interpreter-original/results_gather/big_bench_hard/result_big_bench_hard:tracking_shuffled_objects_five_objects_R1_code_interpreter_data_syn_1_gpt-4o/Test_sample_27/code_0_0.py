# Initial positions
positions = {
    "Alice": "benchwarmer",
    "Bob": "right midfielder",
    "Claire": "cheerleader",
    "Dave": "striker",
    "Eve": "left winger"
}

# Swaps
swaps = [
    ("Alice", "Eve"),
    ("Claire", "Alice"),
    ("Eve", "Alice"),
    ("Bob", "Dave"),
    ("Dave", "Eve")
]

# Apply swaps
for player1, player2 in swaps:
    positions[player1], positions[player2] = positions[player2], positions[player1]

# Output Claire's final position
print(positions["Claire"])