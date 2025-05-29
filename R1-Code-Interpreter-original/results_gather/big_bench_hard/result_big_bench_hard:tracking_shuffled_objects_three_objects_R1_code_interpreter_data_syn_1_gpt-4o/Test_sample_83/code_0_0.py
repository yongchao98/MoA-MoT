# Initial positions
positions = {
    "Alice": "Center Midfielder",
    "Bob": "Left Midfielder",
    "Claire": "Left Winger"
}

# First swap: Alice and Claire
positions["Alice"], positions["Claire"] = positions["Claire"], positions["Alice"]

# Second swap: Claire and Bob
positions["Claire"], positions["Bob"] = positions["Bob"], positions["Claire"]

# Third swap: Bob and Alice
positions["Bob"], positions["Alice"] = positions["Alice"], positions["Bob"]

# Claire's final position
print(positions["Claire"])