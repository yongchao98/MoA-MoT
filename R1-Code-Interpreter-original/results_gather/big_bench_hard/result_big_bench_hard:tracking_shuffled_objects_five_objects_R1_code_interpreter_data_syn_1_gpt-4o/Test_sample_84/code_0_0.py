# Initial book ownership
books = {
    "Alice": "Hound of the Baskervilles",
    "Bob": "Moby Dick",
    "Claire": "The Fellowship of the Ring",
    "Dave": "The Odyssey",
    "Eve": "Frankenstein"
}

# Swap 1: Eve and Alice
books["Alice"], books["Eve"] = books["Eve"], books["Alice"]

# Swap 2: Bob and Eve
books["Bob"], books["Eve"] = books["Eve"], books["Bob"]

# Swap 3: Claire and Dave
books["Claire"], books["Dave"] = books["Dave"], books["Claire"]

# Swap 4: Alice and Eve
books["Alice"], books["Eve"] = books["Eve"], books["Alice"]

# Swap 5: Claire and Bob
books["Claire"], books["Bob"] = books["Bob"], books["Claire"]

# Output the final book Claire has
print(books["Claire"])