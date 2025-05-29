# Initial book assignment
books = {
    "Alice": "Catch-22",
    "Bob": "The Fellowship of the Ring",
    "Claire": "Lolita",
    "Dave": "The Great Gatsby",
    "Eve": "Frankenstein",
    "Fred": "The Odyssey",
    "Gertrude": "Moby Dick"
}

# List of swaps
swaps = [
    ("Fred", "Alice"),
    ("Fred", "Eve"),
    ("Fred", "Dave"),
    ("Claire", "Fred"),
    ("Gertrude", "Dave"),
    ("Fred", "Alice"),
    ("Fred", "Bob")
]

# Perform the swaps
for person1, person2 in swaps:
    books[person1], books[person2] = books[person2], books[person1]

# Output the book Bob has at the end
print(books["Bob"])