# Initial gifts
alice = "yellow"
bob = "brown"
claire = "blue"

# First swap (Bob and Alice)
alice, bob = bob, alice

# Second swap (Claire and Alice)
alice, claire = claire, alice

# Third swap (Bob and Alice)
alice, bob = bob, alice

# Output the final gift Claire has
print(claire)