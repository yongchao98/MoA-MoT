# Initial gifts
alice = "black ball"
bob = "orange ball"
claire = "white present"

# First swap: Claire and Bob
bob, claire = claire, bob

# Second swap: Bob and Alice
alice, bob = bob, alice

# Third swap: Bob and Claire
bob, claire = claire, bob

# Result
print(claire)