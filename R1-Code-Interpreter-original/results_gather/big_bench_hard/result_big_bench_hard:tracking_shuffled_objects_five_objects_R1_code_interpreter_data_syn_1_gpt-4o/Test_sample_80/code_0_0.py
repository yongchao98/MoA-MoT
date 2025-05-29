# Initial positions
positions = {
    "Alice": "right midfielder",
    "Bob": "benchwarmer",
    "Claire": "left midfielder",
    "Dave": "right winger",
    "Eve": "striker"
}

# Swaps
swaps = [
    ("Alice", "Bob"),
    ("Eve", "Claire"),
    ("Eve", "Bob"),
    ("Dave", "Bob"),
    ("Claire", "Alice")
]

# Apply swaps
for player1, player2 in swaps:
    positions[player1], positions[player2] = positions[player2], positions[player1]

# Output Bob's final position
print(positions["Bob"])