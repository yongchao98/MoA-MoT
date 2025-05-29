# Initial positions
positions = {
    "Alice": "goalkeeper",
    "Bob": "center midfielder",
    "Claire": "benchwarmer",
    "Dave": "left winger",
    "Eve": "striker",
    "Fred": "right midfielder",
    "Gertrude": "left midfielder"
}

# Swaps
swaps = [
    ("Eve", "Claire"),
    ("Claire", "Dave"),
    ("Dave", "Bob"),
    ("Fred", "Gertrude"),
    ("Bob", "Claire"),
    ("Dave", "Alice"),
    ("Claire", "Eve")
]

# Apply swaps
for player1, player2 in swaps:
    positions[player1], positions[player2] = positions[player2], positions[player1]

# Determine Dave's final position
dave_final_position = positions["Dave"]
print(dave_final_position)