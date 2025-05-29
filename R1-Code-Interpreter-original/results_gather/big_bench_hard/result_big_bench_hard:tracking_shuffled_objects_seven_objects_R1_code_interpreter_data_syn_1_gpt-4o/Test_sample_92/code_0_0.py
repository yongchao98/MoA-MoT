# Initial positions
positions = {
    "Alice": "fullback",
    "Bob": "benchwarmer",
    "Claire": "right midfielder",
    "Dave": "striker",
    "Eve": "left winger",
    "Fred": "goalkeeper",
    "Gertrude": "cheerleader"
}

# List of swaps
swaps = [
    ("Eve", "Bob"),
    ("Fred", "Claire"),
    ("Bob", "Fred"),
    ("Bob", "Claire"),
    ("Eve", "Alice"),
    ("Dave", "Gertrude"),
    ("Claire", "Gertrude")
]

# Apply swaps
for player1, player2 in swaps:
    positions[player1], positions[player2] = positions[player2], positions[player1]

# Find Dave's final position
dave_final_position = positions["Dave"]
print(dave_final_position)