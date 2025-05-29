# Initial gifts
alice_gift = "yellow present"
bob_gift = "orange ball"
claire_gift = "green present"

# First swap: Claire and Bob
claire_gift, bob_gift = bob_gift, claire_gift

# Second swap: Alice and Claire
alice_gift, claire_gift = claire_gift, alice_gift

# Third swap: Bob and Alice
bob_gift, alice_gift = alice_gift, bob_gift

# Final gift with Alice
print(alice_gift)