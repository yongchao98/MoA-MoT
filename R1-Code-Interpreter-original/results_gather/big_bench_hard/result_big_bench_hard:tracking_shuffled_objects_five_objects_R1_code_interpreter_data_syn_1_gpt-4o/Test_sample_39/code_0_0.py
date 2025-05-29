# Initial book assignments
books = {
    "Alice": "The Great Gatsby",
    "Bob": "The Pearl",
    "Claire": "Ulysses",
    "Dave": "Hound of the Baskervilles",
    "Eve": "Moby Dick"
}

# Sequence of swaps
swaps = [
    ("Dave", "Claire"),
    ("Claire", "Eve"),
    ("Alice", "Claire"),
    ("Claire", "Bob"),
    ("Dave", "Eve")
]

# Perform the swaps
for person1, person2 in swaps:
    books[person1], books[person2] = books[person2], books[person1]

# Output the book Claire has at the end
print(books["Claire"])