# Initial book ownership
alice_book = "The Odyssey"
bob_book = "The Fellowship of the Ring"
claire_book = "Hound of the Baskervilles"

# First swap (Bob and Alice)
alice_book, bob_book = bob_book, alice_book

# Second swap (Claire and Bob)
claire_book, bob_book = bob_book, claire_book

# Third swap (Alice and Claire)
alice_book, claire_book = claire_book, alice_book

# Final book with Bob
print(bob_book)