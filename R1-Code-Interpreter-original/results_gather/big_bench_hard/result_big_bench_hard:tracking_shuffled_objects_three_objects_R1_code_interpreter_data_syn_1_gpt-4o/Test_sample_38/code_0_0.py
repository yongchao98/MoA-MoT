# Initial state
alice = "pink"
bob = "brown"
claire = "orange"

# First swap: Claire and Alice
alice, claire = claire, alice

# Second swap: Alice and Bob
alice, bob = bob, alice

# Third swap: Bob and Claire
bob, claire = claire, bob

# Output the final ball Alice has
print(alice)