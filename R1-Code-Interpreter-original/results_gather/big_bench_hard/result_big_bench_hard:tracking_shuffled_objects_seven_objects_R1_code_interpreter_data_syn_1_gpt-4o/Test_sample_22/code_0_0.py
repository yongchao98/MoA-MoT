# Initial positions
positions = {
    "Alice": "left winger",
    "Bob": "cheerleader",
    "Claire": "fullback",
    "Dave": "right winger",
    "Eve": "striker",
    "Fred": "right midfielder",
    "Gertrude": "center midfielder"
}

# Swaps
swaps = [
    ("Bob", "Dave"),
    ("Eve", "Claire"),
    ("Gertrude", "Alice"),
    ("Bob", "Claire"),
    ("Alice", "Bob"),
    ("Eve", "Fred"),
    ("Alice", "Fred")
]

# Apply swaps
for player1, player2 in swaps:
    positions[player1], positions[player2] = positions[player2], positions[player1]

# Output Eve's final position
print(positions["Eve"])