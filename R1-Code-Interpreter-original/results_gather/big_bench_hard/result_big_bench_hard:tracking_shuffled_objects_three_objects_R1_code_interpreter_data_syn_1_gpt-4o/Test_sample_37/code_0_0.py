# Initial gifts
alice = "orange ball"
bob = "purple present"
claire = "blue present"

# First swap: Claire and Bob
claire, bob = bob, claire

# Second swap: Claire and Alice
claire, alice = alice, claire

# Third swap: Bob and Alice
bob, alice = alice, bob

# Final gift with Alice
print(alice)