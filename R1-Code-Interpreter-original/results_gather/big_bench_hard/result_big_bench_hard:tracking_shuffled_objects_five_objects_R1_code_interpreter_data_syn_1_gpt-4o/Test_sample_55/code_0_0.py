# Initial positions
positions = {
    "Alice": "Right Winger",
    "Bob": "Left Midfielder",
    "Claire": "Right Midfielder",
    "Dave": "Striker",
    "Eve": "Fullback"
}

# Swaps
swaps = [
    ("Eve", "Claire"),
    ("Dave", "Bob"),
    ("Eve", "Dave"),
    ("Alice", "Claire"),
    ("Bob", "Eve")
]

# Apply swaps
for player1, player2 in swaps:
    positions[player1], positions[player2] = positions[player2], positions[player1]

# Output Eve's final position
print(positions["Eve"])