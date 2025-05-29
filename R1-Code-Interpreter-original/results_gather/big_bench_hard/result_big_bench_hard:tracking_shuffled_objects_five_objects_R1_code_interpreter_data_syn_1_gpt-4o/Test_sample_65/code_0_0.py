# Initial book assignment
books = {
    "Alice": "The Great Gatsby",
    "Bob": "Moby Dick",
    "Claire": "The Pearl",
    "Dave": "Lolita",
    "Eve": "The Odyssey"
}

# List of swaps
swaps = [
    ("Eve", "Claire"),
    ("Dave", "Alice"),
    ("Dave", "Eve"),
    ("Alice", "Dave"),
    ("Bob", "Eve")
]

# Perform the swaps
for person1, person2 in swaps:
    books[person1], books[person2] = books[person2], books[person1]

# Output the final book Bob has
print(books["Bob"])