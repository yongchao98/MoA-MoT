# Initial book ownership
books = {
    "Alice": "Ulysses",
    "Bob": "Moby Dick",
    "Claire": "The Pearl",
    "Dave": "The Odyssey",
    "Eve": "The Fellowship of the Ring"
}

# Perform the swaps
# Alice and Claire swap
books["Alice"], books["Claire"] = books["Claire"], books["Alice"]

# Eve and Claire swap
books["Eve"], books["Claire"] = books["Claire"], books["Eve"]

# Claire and Bob swap
books["Claire"], books["Bob"] = books["Bob"], books["Claire"]

# Dave and Claire swap
books["Dave"], books["Claire"] = books["Claire"], books["Dave"]

# Bob and Dave swap
books["Bob"], books["Dave"] = books["Dave"], books["Bob"]

# Output the book Dave has at the end
print(books["Dave"])