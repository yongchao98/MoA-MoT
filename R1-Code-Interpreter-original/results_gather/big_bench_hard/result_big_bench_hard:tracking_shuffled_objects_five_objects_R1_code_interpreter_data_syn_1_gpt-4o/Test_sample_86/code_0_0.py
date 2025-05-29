# Initial book assignment
books = {
    "Alice": "The Great Gatsby",
    "Bob": "The Pearl",
    "Claire": "Hound of the Baskervilles",
    "Dave": "The Fellowship of the Ring",
    "Eve": "The Odyssey"
}

# Swaps
# Claire and Eve swap books
books["Claire"], books["Eve"] = books["Eve"], books["Claire"]

# Alice and Dave swap books
books["Alice"], books["Dave"] = books["Dave"], books["Alice"]

# Eve and Alice swap books
books["Eve"], books["Alice"] = books["Alice"], books["Eve"]

# Dave and Bob swap books
books["Dave"], books["Bob"] = books["Bob"], books["Dave"]

# Dave and Alice swap books
books["Dave"], books["Alice"] = books["Alice"], books["Dave"]

# Output the book Eve has at the end
print(books["Eve"])