# Initial book ownership
books = {
    "Alice": "Ulysses",
    "Bob": "The Odyssey",
    "Claire": "Hound of the Baskervilles",
    "Dave": "Moby Dick",
    "Eve": "Frankenstein",
    "Fred": "The Pearl",
    "Gertrude": "Lolita"
}

# List of swaps
swaps = [
    ("Alice", "Gertrude"),
    ("Dave", "Fred"),
    ("Alice", "Claire"),
    ("Claire", "Bob"),
    ("Eve", "Claire"),
    ("Fred", "Alice"),
    ("Claire", "Dave")
]

# Perform the swaps
for person1, person2 in swaps:
    books[person1], books[person2] = books[person2], books[person1]

# Output the final book with Alice
print(books["Alice"])