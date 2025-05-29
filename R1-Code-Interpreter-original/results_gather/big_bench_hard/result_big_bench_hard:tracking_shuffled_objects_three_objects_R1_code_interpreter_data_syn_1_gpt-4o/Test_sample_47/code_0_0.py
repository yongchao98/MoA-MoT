# Initial book ownership
alice_book = "Frankenstein"
bob_book = "The Odyssey"
claire_book = "The Great Gatsby"

# First swap: Alice and Claire
alice_book, claire_book = claire_book, alice_book

# Second swap: Alice and Bob
alice_book, bob_book = bob_book, alice_book

# Third swap: Alice and Claire
alice_book, claire_book = claire_book, alice_book

# Final book with Claire
print(claire_book)