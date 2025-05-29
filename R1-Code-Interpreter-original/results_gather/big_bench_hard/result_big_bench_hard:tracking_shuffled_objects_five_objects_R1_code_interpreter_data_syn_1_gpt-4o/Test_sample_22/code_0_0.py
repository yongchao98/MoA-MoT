# Initial book ownership
books = {
    "Alice": "Frankenstein",
    "Bob": "Ulysses",
    "Claire": "Lolita",
    "Dave": "The Odyssey",
    "Eve": "The Fellowship of the Ring"
}

# Trades
# Eve and Dave swap books
books["Eve"], books["Dave"] = books["Dave"], books["Eve"]

# Claire and Alice swap books
books["Claire"], books["Alice"] = books["Alice"], books["Claire"]

# Dave and Eve swap books again
books["Dave"], books["Eve"] = books["Eve"], books["Dave"]

# Bob and Claire swap books
books["Bob"], books["Claire"] = books["Claire"], books["Bob"]

# Alice and Claire swap books
books["Alice"], books["Claire"] = books["Claire"], books["Alice"]

# Output the book Eve has at the end
print(books["Eve"])