# Initial book ownership
books = {
    "Alice": "Moby Dick",
    "Bob": "Lolita",
    "Claire": "The Great Gatsby",
    "Dave": "Catch-22",
    "Eve": "Ulysses"
}

# Swap operations
swaps = [
    ("Claire", "Dave"),
    ("Bob", "Claire"),
    ("Dave", "Alice"),
    ("Dave", "Claire"),
    ("Claire", "Eve")
]

# Perform swaps
for person1, person2 in swaps:
    books[person1], books[person2] = books[person2], books[person1]

# Output the book Eve has at the end
print(books["Eve"])