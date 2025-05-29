# Initial ball distribution
balls = {
    "Alice": "green",
    "Bob": "brown",
    "Claire": "orange",
    "Dave": "red",
    "Eve": "blue",
    "Fred": "black",
    "Gertrude": "yellow"
}

# Swaps
swaps = [
    ("Eve", "Dave"),
    ("Gertrude", "Eve"),
    ("Bob", "Fred"),
    ("Claire", "Fred"),
    ("Alice", "Claire"),
    ("Eve", "Fred"),
    ("Eve", "Dave")
]

# Perform swaps
for player1, player2 in swaps:
    balls[player1], balls[player2] = balls[player2], balls[player1]

# Output the final ball Dave has
print(balls["Dave"])