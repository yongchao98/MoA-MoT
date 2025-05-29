# Initial book ownership
books = {
    "Alice": "The Fellowship of the Ring",
    "Bob": "The Odyssey",
    "Claire": "Frankenstein",
    "Dave": "Hound of the Baskervilles",
    "Eve": "Ulysses"
}

# Perform the swaps
# 1. Alice and Claire swap
books["Alice"], books["Claire"] = books["Claire"], books["Alice"]

# 2. Alice and Eve swap
books["Alice"], books["Eve"] = books["Eve"], books["Alice"]

# 3. Dave and Claire swap
books["Dave"], books["Claire"] = books["Claire"], books["Dave"]

# 4. Dave and Bob swap
books["Dave"], books["Bob"] = books["Bob"], books["Dave"]

# 5. Dave and Alice swap
books["Dave"], books["Alice"] = books["Alice"], books["Dave"]

# Output the book Alice has at the end
print(books["Alice"])