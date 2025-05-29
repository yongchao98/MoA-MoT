# Initial book ownership
books = {
    "Alice": "The Great Gatsby",
    "Bob": "Frankenstein",
    "Claire": "Moby Dick",
    "Dave": "The Fellowship of the Ring",
    "Eve": "Catch-22"
}

# Trades
# 1. Alice and Bob swap books
books["Alice"], books["Bob"] = books["Bob"], books["Alice"]

# 2. Dave and Eve swap books
books["Dave"], books["Eve"] = books["Eve"], books["Dave"]

# 3. Eve and Alice swap books
books["Eve"], books["Alice"] = books["Alice"], books["Eve"]

# 4. Dave and Eve swap books
books["Dave"], books["Eve"] = books["Eve"], books["Dave"]

# 5. Claire and Dave swap books
books["Claire"], books["Dave"] = books["Dave"], books["Claire"]

# Output the book Eve has at the end
print(books["Eve"])