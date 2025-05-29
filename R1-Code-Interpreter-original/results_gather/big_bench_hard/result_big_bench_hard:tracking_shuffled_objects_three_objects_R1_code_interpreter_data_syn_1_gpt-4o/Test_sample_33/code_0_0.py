# Initial gifts
alice = "white present"
bob = "pink ball"
claire = "black ball"

# First swap: Alice and Claire
alice, claire = claire, alice

# Second swap: Alice and Bob
alice, bob = bob, alice

# Third swap: Claire and Alice
claire, alice = alice, claire

# Output what Claire has at the end
print(claire)