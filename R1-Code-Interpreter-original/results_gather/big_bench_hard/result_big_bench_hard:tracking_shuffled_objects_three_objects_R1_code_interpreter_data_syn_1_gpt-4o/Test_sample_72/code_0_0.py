# Initial ball possession
alice = "yellow"
bob = "purple"
claire = "blue"

# First swap between Alice and Claire
alice, claire = claire, alice

# Second swap between Bob and Claire
bob, claire = claire, bob

# Third swap between Alice and Claire
alice, claire = claire, alice

# Output the final ball with Bob
print(bob)