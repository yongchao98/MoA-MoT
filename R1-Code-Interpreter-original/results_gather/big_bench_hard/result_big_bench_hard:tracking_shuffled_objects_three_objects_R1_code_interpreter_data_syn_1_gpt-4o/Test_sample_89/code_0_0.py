# Initial gifts
alice = "black ball"
bob = "green present"
claire = "brown present"

# First swap: Alice and Bob
alice, bob = bob, alice

# Second swap: Claire and Bob
claire, bob = bob, claire

# Third swap: Alice and Claire
alice, claire = claire, alice

# Output the final gift Alice has
print(alice)