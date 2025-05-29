# Initial book assignment
books = {
    "Alice": "Hound of the Baskervilles",
    "Bob": "Moby Dick",
    "Claire": "The Fellowship of the Ring",
    "Dave": "Catch-22",
    "Eve": "Frankenstein"
}

# Perform the swaps
# Claire and Dave swap
books["Claire"], books["Dave"] = books["Dave"], books["Claire"]

# Alice and Eve swap
books["Alice"], books["Eve"] = books["Eve"], books["Alice"]

# Bob and Eve swap
books["Bob"], books["Eve"] = books["Eve"], books["Bob"]

# Dave and Bob swap
books["Dave"], books["Bob"] = books["Bob"], books["Dave"]

# Alice and Dave swap
books["Alice"], books["Dave"] = books["Dave"], books["Alice"]

# Output the book Dave has at the end
print(books["Dave"])