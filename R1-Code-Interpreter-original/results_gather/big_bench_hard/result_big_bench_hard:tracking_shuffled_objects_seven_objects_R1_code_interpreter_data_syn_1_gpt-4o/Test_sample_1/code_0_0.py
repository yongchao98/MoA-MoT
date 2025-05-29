# Initial ball assignment
balls = {
    "Alice": "orange",
    "Bob": "yellow",
    "Claire": "brown",
    "Dave": "white",
    "Eve": "black",
    "Fred": "red",
    "Gertrude": "purple"
}

# Swaps
swaps = [
    ("Claire", "Bob"),
    ("Claire", "Dave"),
    ("Fred", "Gertrude"),
    ("Eve", "Claire"),
    ("Eve", "Fred"),
    ("Alice", "Dave"),
    ("Eve", "Claire")
]

# Perform swaps
for player1, player2 in swaps:
    balls[player1], balls[player2] = balls[player2], balls[player1]

# Output the ball Eve has at the end
print(balls["Eve"])