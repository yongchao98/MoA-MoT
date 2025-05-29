# Initial book assignments
books = {
    "Alice": "The Great Gatsby",
    "Bob": "The Fellowship of the Ring",
    "Claire": "Ulysses",
    "Dave": "The Odyssey",
    "Eve": "Frankenstein",
    "Fred": "Hound of the Baskervilles",
    "Gertrude": "Moby Dick"
}

# List of trades
trades = [
    ("Claire", "Eve"),
    ("Dave", "Fred"),
    ("Alice", "Eve"),
    ("Claire", "Dave"),
    ("Fred", "Claire"),
    ("Claire", "Gertrude"),
    ("Fred", "Bob")
]

# Perform the trades
for person1, person2 in trades:
    books[person1], books[person2] = books[person2], books[person1]

# Output the book Alice has at the end
print(books["Alice"])