# Initial gifts
alice = "pink ball"
bob = "yellow present"
claire = "orange ball"

# First swap (Bob and Claire)
bob, claire = claire, bob

# Second swap (Alice and Claire)
alice, claire = claire, alice

# Third swap (Bob and Claire)
bob, claire = claire, bob

# Output the final gift Bob has
print(bob)