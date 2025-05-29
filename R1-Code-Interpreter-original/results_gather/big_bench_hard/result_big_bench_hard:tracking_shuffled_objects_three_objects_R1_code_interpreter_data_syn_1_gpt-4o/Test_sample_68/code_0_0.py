# Initial gifts
alice = "pink ball"
bob = "blue present"
claire = "brown present"

# First swap: Claire and Alice
alice, claire = claire, alice

# Second swap: Bob and Alice
alice, bob = bob, alice

# Third swap: Alice and Claire
alice, claire = claire, alice

# Final gift with Alice
print(alice)