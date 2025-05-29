# Initial distribution of presents
presents = {
    "Alice": "blue",
    "Bob": "yellow",
    "Claire": "red",
    "Dave": "black",
    "Eve": "white",
    "Fred": "brown",
    "Gertrude": "orange"
}

# List of swaps
swaps = [
    ("Alice", "Fred"),
    ("Claire", "Bob"),
    ("Dave", "Fred"),
    ("Eve", "Alice"),
    ("Bob", "Alice"),
    ("Eve", "Gertrude"),
    ("Fred", "Alice")
]

# Perform the swaps
for person1, person2 in swaps:
    presents[person1], presents[person2] = presents[person2], presents[person1]

# Output the present Bob has at the end
print(presents["Bob"])