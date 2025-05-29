# Initial book assignment
books = {
    "Alice": "Moby Dick",
    "Bob": "Lolita",
    "Claire": "Ulysses",
    "Dave": "The Fellowship of the Ring",
    "Eve": "Hound of the Baskervilles",
    "Fred": "Catch-22",
    "Gertrude": "Frankenstein"
}

# Perform the swaps
swaps = [
    ("Fred", "Eve"),
    ("Gertrude", "Alice"),
    ("Dave", "Claire"),
    ("Alice", "Bob"),
    ("Alice", "Fred"),
    ("Eve", "Dave"),
    ("Eve", "Alice")
]

for person1, person2 in swaps:
    books[person1], books[person2] = books[person2], books[person1]

# Output the book Gertrude has at the end
print(books["Gertrude"])