# Initial book ownership
alice_book = "Hound of the Baskervilles"
bob_book = "The Odyssey"
claire_book = "Catch-22"

# First swap: Claire and Alice
alice_book, claire_book = claire_book, alice_book

# Second swap: Bob and Claire
bob_book, claire_book = claire_book, bob_book

# Third swap: Alice and Bob
alice_book, bob_book = bob_book, alice_book

# Final book with Alice
print(alice_book)