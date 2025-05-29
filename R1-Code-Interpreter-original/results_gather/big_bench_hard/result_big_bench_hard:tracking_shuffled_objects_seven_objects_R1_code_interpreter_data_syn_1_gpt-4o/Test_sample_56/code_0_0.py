# Initial distribution of presents
presents = {
    "Alice": "purple",
    "Bob": "orange",
    "Claire": "white",
    "Dave": "green",
    "Eve": "yellow",
    "Fred": "brown",
    "Gertrude": "red"
}

# List of swaps
swaps = [
    ("Eve", "Bob"),
    ("Dave", "Claire"),
    ("Alice", "Bob"),
    ("Alice", "Gertrude"),
    ("Claire", "Bob"),
    ("Dave", "Fred"),
    ("Bob", "Eve")
]

# Perform the swaps
for person1, person2 in swaps:
    presents[person1], presents[person2] = presents[person2], presents[person1]

# Output the present Claire has at the end
print(presents["Claire"])