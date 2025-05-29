# Initial positions
positions = {
    "Alice": "left winger",
    "Bob": "goalkeeper",
    "Claire": "right winger",
    "Dave": "fullback",
    "Eve": "right midfielder"
}

# Swaps
swaps = [
    ("Claire", "Alice"),
    ("Claire", "Dave"),
    ("Bob", "Claire"),
    ("Eve", "Alice"),
    ("Claire", "Bob")
]

# Perform swaps
for player1, player2 in swaps:
    positions[player1], positions[player2] = positions[player2], positions[player1]

# Output Eve's final position
print(positions["Eve"])