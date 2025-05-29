# Initial state
alice = "yellow"
bob = "purple"
claire = "blue"

# First swap: Alice and Claire
alice, claire = claire, alice  # Alice gets blue, Claire gets yellow

# Second swap: Bob and Claire
bob, claire = claire, bob  # Bob gets yellow, Claire gets purple

# Third swap: Alice and Claire
alice, claire = claire, alice  # Alice gets purple, Claire gets blue

# Output the final ball Alice has
print(alice)