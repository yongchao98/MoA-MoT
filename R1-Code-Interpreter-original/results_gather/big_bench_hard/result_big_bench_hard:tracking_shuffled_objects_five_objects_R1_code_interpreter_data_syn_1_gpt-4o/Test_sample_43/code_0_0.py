# Initial book distribution
books = {
    "Alice": "The Fellowship of the Ring",
    "Bob": "The Odyssey",
    "Claire": "Catch-22",
    "Dave": "The Great Gatsby",
    "Eve": "Lolita"
}

# Perform the trades
# 1. Alice and Eve swap books
books["Alice"], books["Eve"] = books["Eve"], books["Alice"]

# 2. Claire and Eve swap books
books["Claire"], books["Eve"] = books["Eve"], books["Claire"]

# 3. Alice and Eve swap books
books["Alice"], books["Eve"] = books["Eve"], books["Alice"]

# 4. Eve and Dave swap books
books["Eve"], books["Dave"] = books["Dave"], books["Eve"]

# 5. Alice and Bob swap books
books["Alice"], books["Bob"] = books["Bob"], books["Alice"]

# Output the book Bob has at the end
print(books["Bob"])