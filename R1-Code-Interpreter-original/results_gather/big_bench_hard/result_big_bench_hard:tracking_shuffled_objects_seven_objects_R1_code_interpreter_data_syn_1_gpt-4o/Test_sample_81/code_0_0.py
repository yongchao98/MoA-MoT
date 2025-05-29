# Initial positions
positions = {
    "Alice": "right winger",
    "Bob": "right midfielder",
    "Claire": "cheerleader",
    "Dave": "striker",
    "Eve": "center midfielder",
    "Fred": "benchwarmer",
    "Gertrude": "left midfielder"
}

# Swaps
swaps = [
    ("Eve", "Claire"),
    ("Gertrude", "Dave"),
    ("Bob", "Claire"),
    ("Fred", "Dave"),
    ("Gertrude", "Fred"),
    ("Fred", "Claire"),
    ("Claire", "Alice")
]

# Apply swaps
for player1, player2 in swaps:
    positions[player1], positions[player2] = positions[player2], positions[player1]

# Output Claire's final position
print(positions["Claire"])