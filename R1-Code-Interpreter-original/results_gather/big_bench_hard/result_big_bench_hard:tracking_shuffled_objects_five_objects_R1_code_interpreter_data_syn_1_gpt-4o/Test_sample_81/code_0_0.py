# Initial book assignments
books = {
    "Alice": "The Great Gatsby",
    "Bob": "The Odyssey",
    "Claire": "The Fellowship of the Ring",
    "Dave": "Moby Dick",
    "Eve": "Lolita"
}

# Perform the trades
# Trade 1: Alice and Bob swap books
books["Alice"], books["Bob"] = books["Bob"], books["Alice"]

# Trade 2: Eve and Bob swap books
books["Eve"], books["Bob"] = books["Bob"], books["Eve"]

# Trade 3: Bob and Claire swap books
books["Bob"], books["Claire"] = books["Claire"], books["Bob"]

# Trade 4: Alice and Claire swap books
books["Alice"], books["Claire"] = books["Claire"], books["Alice"]

# Trade 5: Dave and Claire swap books
books["Dave"], books["Claire"] = books["Claire"], books["Dave"]

# Output the book Dave has at the end
print(books["Dave"])