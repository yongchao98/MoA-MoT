# Initial book distribution
books = {
    "Alice": "Frankenstein",
    "Bob": "The Pearl",
    "Claire": "Hound of the Baskervilles",
    "Dave": "Ulysses",
    "Eve": "The Great Gatsby"
}

# Swap 1: Dave and Claire
books["Dave"], books["Claire"] = books["Claire"], books["Dave"]

# Swap 2: Bob and Claire
books["Bob"], books["Claire"] = books["Claire"], books["Bob"]

# Swap 3: Claire and Eve
books["Claire"], books["Eve"] = books["Eve"], books["Claire"]

# Swap 4: Bob and Claire
books["Bob"], books["Claire"] = books["Claire"], books["Bob"]

# Swap 5: Alice and Claire
books["Alice"], books["Claire"] = books["Claire"], books["Alice"]

# Output the book Dave has at the end
print(books["Dave"])