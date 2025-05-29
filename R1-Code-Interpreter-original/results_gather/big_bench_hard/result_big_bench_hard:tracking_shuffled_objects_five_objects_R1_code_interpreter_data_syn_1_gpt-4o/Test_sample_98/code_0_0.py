# Initial book ownership
books = {
    "Alice": "The Great Gatsby",
    "Bob": "Frankenstein",
    "Claire": "The Pearl",
    "Dave": "Catch-22",
    "Eve": "Ulysses"
}

# Sequence of swaps
swaps = [
    ("Claire", "Dave"),
    ("Bob", "Alice"),
    ("Dave", "Alice"),
    ("Bob", "Claire"),
    ("Dave", "Eve")
]

# Perform the swaps
for person1, person2 in swaps:
    books[person1], books[person2] = books[person2], books[person1]

# Output the book Eve has at the end
print(books["Eve"])