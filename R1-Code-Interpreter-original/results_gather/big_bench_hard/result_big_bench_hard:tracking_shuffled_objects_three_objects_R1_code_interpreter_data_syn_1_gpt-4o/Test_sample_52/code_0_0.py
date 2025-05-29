# Initial state
alice = "Brown"
bob = "Blue"
claire = "Purple"

# First swap: Bob and Alice
alice, bob = bob, alice

# Second swap: Claire and Alice
alice, claire = claire, alice

# Third swap: Bob and Claire
bob, claire = claire, bob

# Final state of Alice
print(alice)