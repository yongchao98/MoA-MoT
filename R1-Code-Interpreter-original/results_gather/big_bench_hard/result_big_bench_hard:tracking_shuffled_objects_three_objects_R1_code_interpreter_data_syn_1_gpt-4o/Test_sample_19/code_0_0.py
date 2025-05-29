# Initial book ownership
alice_book = "The Great Gatsby"
bob_book = "Frankenstein"
claire_book = "Ulysses"

# First swap: Alice and Claire
alice_book, claire_book = claire_book, alice_book

# Second swap: Bob and Claire
bob_book, claire_book = claire_book, bob_book

# Third swap: Claire and Alice
claire_book, alice_book = alice_book, claire_book

# Final book with Alice
print(alice_book)