# Initial book assignment
books = {
    "Alice": "Frankenstein",
    "Bob": "Ulysses",
    "Claire": "The Great Gatsby",
    "Dave": "The Odyssey",
    "Eve": "The Fellowship of the Ring",
    "Fred": "Catch-22",
    "Gertrude": "Moby Dick"
}

# Perform the trades
# Claire and Gertrude swap books
books["Claire"], books["Gertrude"] = books["Gertrude"], books["Claire"]

# Fred and Dave swap books
books["Fred"], books["Dave"] = books["Dave"], books["Fred"]

# Bob and Alice swap books
books["Bob"], books["Alice"] = books["Alice"], books["Bob"]

# Eve and Gertrude swap books
books["Eve"], books["Gertrude"] = books["Gertrude"], books["Eve"]

# Dave and Bob swap books
books["Dave"], books["Bob"] = books["Bob"], books["Dave"]

# Fred and Eve swap books
books["Fred"], books["Eve"] = books["Eve"], books["Fred"]

# Bob and Fred swap books
books["Bob"], books["Fred"] = books["Fred"], books["Bob"]

# Output the final book Alice has
print(books["Alice"])