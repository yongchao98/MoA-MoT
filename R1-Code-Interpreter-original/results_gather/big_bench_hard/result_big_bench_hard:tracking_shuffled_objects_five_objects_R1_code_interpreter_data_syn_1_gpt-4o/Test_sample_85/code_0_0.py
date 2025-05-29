# Initial positions
positions = {
    "Alice": "left midfielder",
    "Bob": "goalkeeper",
    "Claire": "fullback",
    "Dave": "striker",
    "Eve": "left winger"
}

# Swaps
swaps = [
    ("Eve", "Dave"),
    ("Dave", "Alice"),
    ("Bob", "Dave"),
    ("Claire", "Bob"),
    ("Claire", "Alice")
]

# Perform swaps
for player1, player2 in swaps:
    positions[player1], positions[player2] = positions[player2], positions[player1]

# Output Bob's final position
print(positions["Bob"])