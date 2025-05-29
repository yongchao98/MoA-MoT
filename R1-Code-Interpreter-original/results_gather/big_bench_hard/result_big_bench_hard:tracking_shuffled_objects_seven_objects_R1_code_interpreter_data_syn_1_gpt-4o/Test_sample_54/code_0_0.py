# Initial book assignments
books = {
    "Alice": "The Odyssey",
    "Bob": "The Great Gatsby",
    "Claire": "Lolita",
    "Dave": "The Pearl",
    "Eve": "The Fellowship of the Ring",
    "Fred": "Frankenstein",
    "Gertrude": "Hound of the Baskervilles"
}

# Perform the swaps
# Claire and Bob swap books
books["Claire"], books["Bob"] = books["Bob"], books["Claire"]

# Alice and Dave swap books
books["Alice"], books["Dave"] = books["Dave"], books["Alice"]

# Gertrude and Fred swap books
books["Gertrude"], books["Fred"] = books["Fred"], books["Gertrude"]

# Bob and Alice swap books
books["Bob"], books["Alice"] = books["Alice"], books["Bob"]

# Alice and Eve swap books
books["Alice"], books["Eve"] = books["Eve"], books["Alice"]

# Dave and Gertrude swap books
books["Dave"], books["Gertrude"] = books["Gertrude"], books["Dave"]

# Dave and Eve swap books
books["Dave"], books["Eve"] = books["Eve"], books["Dave"]

# Output the book Gertrude has at the end
print(books["Gertrude"])