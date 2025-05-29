# Initial book assignment
books = {
    "Alice": "Moby Dick",
    "Bob": "Hound of the Baskervilles",
    "Claire": "The Great Gatsby",
    "Dave": "The Odyssey",
    "Eve": "Ulysses",
    "Fred": "Lolita",
    "Gertrude": "The Pearl"
}

# Swaps
swaps = [
    ("Alice", "Fred"),
    ("Eve", "Fred"),
    ("Fred", "Bob"),
    ("Eve", "Gertrude"),
    ("Gertrude", "Claire"),
    ("Claire", "Bob"),
    ("Dave", "Alice")
]

# Perform swaps
for person1, person2 in swaps:
    books[person1], books[person2] = books[person2], books[person1]

# Output the book Bob has at the end
print(books["Bob"])