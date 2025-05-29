# Initial state
alice = "Pink"
bob = "Brown"
claire = "Orange"

# First swap: Claire and Alice
alice, claire = claire, alice

# Second swap: Alice and Bob
alice, bob = bob, alice

# Third swap: Bob and Claire
bob, claire = claire, bob

# Final state of Claire's ball
print(claire)