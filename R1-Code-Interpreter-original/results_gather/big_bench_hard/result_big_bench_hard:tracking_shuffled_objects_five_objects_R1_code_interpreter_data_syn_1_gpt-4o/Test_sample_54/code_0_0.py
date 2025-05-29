# Initial presents
presents = {
    "Alice": "Purple",
    "Bob": "Blue",
    "Claire": "Brown",
    "Dave": "Yellow",
    "Eve": "Black"
}

# Swap operations
swaps = [
    ("Alice", "Bob"),
    ("Bob", "Dave"),
    ("Claire", "Bob"),
    ("Dave", "Alice"),
    ("Dave", "Eve")
]

# Perform swaps
for person1, person2 in swaps:
    presents[person1], presents[person2] = presents[person2], presents[person1]

# Output the present Dave has at the end
print(presents["Dave"])