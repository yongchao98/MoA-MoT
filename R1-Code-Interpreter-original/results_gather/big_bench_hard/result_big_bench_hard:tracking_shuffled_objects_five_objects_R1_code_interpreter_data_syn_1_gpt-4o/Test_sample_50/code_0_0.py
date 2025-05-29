# Initial book distribution
books = {
    "Alice": "Catch-22",
    "Bob": "The Odyssey",
    "Claire": "Moby Dick",
    "Dave": "Lolita",
    "Eve": "Frankenstein"
}

# Perform the swaps
# Claire and Eve swap books
books["Claire"], books["Eve"] = books["Eve"], books["Claire"]

# Alice and Bob swap books
books["Alice"], books["Bob"] = books["Bob"], books["Alice"]

# Bob and Eve swap books
books["Bob"], books["Eve"] = books["Eve"], books["Bob"]

# Eve and Alice swap books
books["Eve"], books["Alice"] = books["Alice"], books["Eve"]

# Bob and Dave swap books
books["Bob"], books["Dave"] = books["Dave"], books["Bob"]

# Output the book Claire has at the end
print(books["Claire"])