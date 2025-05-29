# Initial state
alice = "pink"
bob = "orange"
claire = "white"

# First swap: Bob and Claire
bob, claire = claire, bob

# Second swap: Claire and Alice
claire, alice = alice, claire

# Third swap: Alice and Bob
alice, bob = bob, alice

# Final ball with Alice
print(alice)