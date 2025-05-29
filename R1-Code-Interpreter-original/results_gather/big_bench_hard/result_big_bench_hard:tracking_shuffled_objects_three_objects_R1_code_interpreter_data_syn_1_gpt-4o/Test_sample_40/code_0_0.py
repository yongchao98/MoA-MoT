# Initial presents
alice_present = "blue"
bob_present = "red"
claire_present = "green"

# First swap: Alice and Bob
alice_present, bob_present = bob_present, alice_present

# Second swap: Claire and Bob
claire_present, bob_present = bob_present, claire_present

# Third swap: Alice and Claire
alice_present, claire_present = claire_present, alice_present

# Result
print(claire_present)