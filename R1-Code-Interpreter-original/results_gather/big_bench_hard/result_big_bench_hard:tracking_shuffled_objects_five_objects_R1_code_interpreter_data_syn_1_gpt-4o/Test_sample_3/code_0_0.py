# Initial book assignment
books = {
    "Alice": "Catch-22",
    "Bob": "Hound of the Baskervilles",
    "Claire": "Frankenstein",
    "Dave": "The Pearl",
    "Eve": "The Fellowship of the Ring"
}

# Perform swaps
# Eve and Alice swap books
books["Alice"], books["Eve"] = books["Eve"], books["Alice"]

# Alice and Claire swap books
books["Alice"], books["Claire"] = books["Claire"], books["Alice"]

# Alice and Bob swap books
books["Alice"], books["Bob"] = books["Bob"], books["Alice"]

# Dave and Alice swap books
books["Dave"], books["Alice"] = books["Alice"], books["Dave"]

# Dave and Claire swap books
books["Dave"], books["Claire"] = books["Claire"], books["Dave"]

# Output the book Dave has at the end
print(books["Dave"])