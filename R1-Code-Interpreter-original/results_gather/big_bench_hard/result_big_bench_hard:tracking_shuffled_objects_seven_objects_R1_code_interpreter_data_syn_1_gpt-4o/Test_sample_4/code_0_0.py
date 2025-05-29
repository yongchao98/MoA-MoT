# Initial book ownership
books = {
    "Alice": "Lolita",
    "Bob": "Frankenstein",
    "Claire": "Moby Dick",
    "Dave": "The Pearl",
    "Eve": "Hound of the Baskervilles",
    "Fred": "The Fellowship of the Ring",
    "Gertrude": "The Odyssey"
}

# Perform the swaps
# 1. Dave and Gertrude swap books
books["Dave"], books["Gertrude"] = books["Gertrude"], books["Dave"]

# 2. Claire and Bob swap books
books["Claire"], books["Bob"] = books["Bob"], books["Claire"]

# 3. Alice and Dave swap books
books["Alice"], books["Dave"] = books["Dave"], books["Alice"]

# 4. Fred and Claire swap books
books["Fred"], books["Claire"] = books["Claire"], books["Fred"]

# 5. Dave and Fred swap books
books["Dave"], books["Fred"] = books["Fred"], books["Dave"]

# 6. Bob and Claire swap books
books["Bob"], books["Claire"] = books["Claire"], books["Bob"]

# 7. Eve and Fred swap books
books["Eve"], books["Fred"] = books["Fred"], books["Eve"]

# Output the book Eve has at the end
print(books["Eve"])