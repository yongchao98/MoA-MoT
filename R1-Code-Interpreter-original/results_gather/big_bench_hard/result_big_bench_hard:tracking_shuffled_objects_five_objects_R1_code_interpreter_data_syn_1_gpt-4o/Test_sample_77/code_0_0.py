# Initial ball positions
balls = {
    "Alice": "black",
    "Bob": "brown",
    "Claire": "blue",
    "Dave": "white",
    "Eve": "green"
}

# Swap sequences
swaps = [
    ("Dave", "Eve"),
    ("Bob", "Claire"),
    ("Claire", "Alice"),
    ("Dave", "Eve"),
    ("Eve", "Claire")
]

# Perform swaps
for p1, p2 in swaps:
    balls[p1], balls[p2] = balls[p2], balls[p1]

# Output the ball Alice has at the end
print(balls["Alice"])