# Initial gifts
alice = "pink ball"
bob = "white present"
claire = "blue present"

# First swap: Bob and Claire
bob, claire = claire, bob

# Second swap: Claire and Alice
claire, alice = alice, claire

# Third swap: Claire and Bob
claire, bob = bob, claire

# Final gift with Alice
print(alice)