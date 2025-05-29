# Initial state
alice = "brown"
bob = "red"
claire = "white"

# First swap: Claire and Bob
claire, bob = bob, claire

# Second swap: Bob and Alice
bob, alice = alice, bob

# Third swap: Claire and Alice
claire, alice = alice, claire

# Final state of Bob
print(bob)