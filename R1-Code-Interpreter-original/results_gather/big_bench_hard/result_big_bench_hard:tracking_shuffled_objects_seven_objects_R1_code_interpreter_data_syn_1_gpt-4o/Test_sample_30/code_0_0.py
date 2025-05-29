# Initial book assignment
books = {
    "Alice": "Ulysses",
    "Bob": "Frankenstein",
    "Claire": "Moby Dick",
    "Dave": "The Odyssey",
    "Eve": "The Great Gatsby",
    "Fred": "Lolita",
    "Gertrude": "The Pearl"
}

# List of swaps
swaps = [
    ("Bob", "Claire"),
    ("Eve", "Alice"),
    ("Alice", "Dave"),
    ("Gertrude", "Dave"),
    ("Alice", "Bob"),
    ("Eve", "Claire"),
    ("Alice", "Fred")
]

# Perform the swaps
for person1, person2 in swaps:
    books[person1], books[person2] = books[person2], books[person1]

# Output the final book Gertrude has
print(books["Gertrude"])