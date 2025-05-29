# Initial gifts
alice = "white present"
bob = "blue present"
claire = "pink ball"

# First swap: Alice and Claire
alice, claire = claire, alice

# Second swap: Claire and Bob
claire, bob = bob, claire

# Third swap: Bob and Alice
bob, alice = alice, bob

# Output the final gift Bob has
print(bob)