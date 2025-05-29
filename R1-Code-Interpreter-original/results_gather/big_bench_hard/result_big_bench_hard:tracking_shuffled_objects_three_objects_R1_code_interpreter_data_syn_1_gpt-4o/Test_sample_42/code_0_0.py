# Initial positions
positions = {
    "Alice": "cheerleader",
    "Bob": "benchwarmer",
    "Claire": "left winger"
}

# First swap: Claire and Bob
positions["Claire"], positions["Bob"] = positions["Bob"], positions["Claire"]

# Second swap: Claire and Alice
positions["Claire"], positions["Alice"] = positions["Alice"], positions["Claire"]

# Third swap: Claire and Bob
positions["Claire"], positions["Bob"] = positions["Bob"], positions["Claire"]

# Final position of Claire
print(positions["Claire"])