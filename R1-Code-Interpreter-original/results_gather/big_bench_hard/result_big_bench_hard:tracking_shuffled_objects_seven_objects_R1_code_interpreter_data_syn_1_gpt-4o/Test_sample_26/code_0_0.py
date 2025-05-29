# Initial book assignment
books = {
    "Alice": "The Great Gatsby",
    "Bob": "Catch-22",
    "Claire": "The Fellowship of the Ring",
    "Dave": "Ulysses",
    "Eve": "Frankenstein",
    "Fred": "The Pearl",
    "Gertrude": "The Odyssey"
}

# Perform the swaps
books["Alice"], books["Fred"] = books["Fred"], books["Alice"]
books["Gertrude"], books["Eve"] = books["Eve"], books["Gertrude"]
books["Alice"], books["Eve"] = books["Eve"], books["Alice"]
books["Dave"], books["Gertrude"] = books["Gertrude"], books["Dave"]
books["Gertrude"], books["Claire"] = books["Claire"], books["Gertrude"]
books["Alice"], books["Claire"] = books["Claire"], books["Alice"]
books["Dave"], books["Bob"] = books["Bob"], books["Dave"]

# Output the final book Dave has
print(books["Dave"])