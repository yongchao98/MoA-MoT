# Initial book assignments
books = {
    "Alice": "The Pearl",
    "Bob": "The Great Gatsby",
    "Claire": "Lolita",
    "Dave": "Ulysses",
    "Eve": "Frankenstein",
    "Fred": "Moby Dick",
    "Gertrude": "The Odyssey"
}

# List of swaps
swaps = [
    ("Gertrude", "Fred"),
    ("Bob", "Dave"),
    ("Eve", "Alice"),
    ("Bob", "Fred"),
    ("Dave", "Eve"),
    ("Claire", "Gertrude"),
    ("Fred", "Eve")
]

# Perform the swaps
for person1, person2 in swaps:
    books[person1], books[person2] = books[person2], books[person1]

# Output the book Dave has at the end
print(books["Dave"])