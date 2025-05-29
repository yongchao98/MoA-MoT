# Initial positions
positions = {
    "Alice": "Left Winger",
    "Bob": "Fullback",
    "Claire": "Right Winger"
}

# First swap: Alice and Claire
positions["Alice"], positions["Claire"] = positions["Claire"], positions["Alice"]

# Second swap: Alice and Bob
positions["Alice"], positions["Bob"] = positions["Bob"], positions["Alice"]

# Third swap: Claire and Alice
positions["Claire"], positions["Alice"] = positions["Alice"], positions["Claire"]

# Output Claire's final position
print(positions["Claire"])