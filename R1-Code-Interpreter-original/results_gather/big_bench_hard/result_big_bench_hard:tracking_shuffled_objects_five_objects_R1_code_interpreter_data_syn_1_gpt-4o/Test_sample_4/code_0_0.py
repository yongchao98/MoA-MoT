# Initial book distribution
books = {
    "Alice": "Moby Dick",
    "Bob": "Ulysses",
    "Claire": "The Fellowship of the Ring",
    "Dave": "Lolita",
    "Eve": "Catch-22"
}

# Swap operations
# Claire and Bob swap books
books["Claire"], books["Bob"] = books["Bob"], books["Claire"]

# Claire and Dave swap books
books["Claire"], books["Dave"] = books["Dave"], books["Claire"]

# Dave and Alice swap books
books["Dave"], books["Alice"] = books["Alice"], books["Dave"]

# Dave and Bob swap books
books["Dave"], books["Bob"] = books["Bob"], books["Dave"]

# Eve and Dave swap books
books["Eve"], books["Dave"] = books["Dave"], books["Eve"]

# Output the book Claire has at the end
print(books["Claire"])