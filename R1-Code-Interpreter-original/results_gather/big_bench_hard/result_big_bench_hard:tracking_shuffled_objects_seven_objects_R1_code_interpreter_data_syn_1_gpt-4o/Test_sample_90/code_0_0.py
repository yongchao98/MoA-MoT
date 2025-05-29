# Initial book assignment
books = {
    "Alice": "The Pearl",
    "Bob": "The Fellowship of the Ring",
    "Claire": "Catch-22",
    "Dave": "The Great Gatsby",
    "Eve": "Frankenstein",
    "Fred": "Ulysses",
    "Gertrude": "The Odyssey"
}

# Perform the swaps
# Gertrude and Bob swap books
books["Gertrude"], books["Bob"] = books["Bob"], books["Gertrude"]

# Alice and Fred swap books
books["Alice"], books["Fred"] = books["Fred"], books["Alice"]

# Gertrude and Claire swap books
books["Gertrude"], books["Claire"] = books["Claire"], books["Gertrude"]

# Fred and Eve swap books
books["Fred"], books["Eve"] = books["Eve"], books["Fred"]

# Claire and Dave swap books
books["Claire"], books["Dave"] = books["Dave"], books["Claire"]

# Eve and Bob swap books
books["Eve"], books["Bob"] = books["Bob"], books["Eve"]

# Fred and Gertrude swap books
books["Fred"], books["Gertrude"] = books["Gertrude"], books["Fred"]

# Output the book Fred has at the end
print(books["Fred"])