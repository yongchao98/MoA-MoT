# Initial distribution of presents
presents = {
    "Alice": "Purple",
    "Bob": "Green",
    "Claire": "Brown",
    "Dave": "Red",
    "Eve": "White"
}

# Sequence of swaps
swaps = [
    ("Eve", "Claire"),
    ("Claire", "Bob"),
    ("Bob", "Alice"),
    ("Eve", "Alice"),
    ("Bob", "Dave")
]

# Perform the swaps
for person1, person2 in swaps:
    presents[person1], presents[person2] = presents[person2], presents[person1]

# Output the present Eve has at the end
print(presents["Eve"])