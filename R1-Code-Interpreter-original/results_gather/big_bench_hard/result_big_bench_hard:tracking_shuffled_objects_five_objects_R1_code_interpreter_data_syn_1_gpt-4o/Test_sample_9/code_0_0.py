# Initial positions
positions = {
    "Alice": "benchwarmer",
    "Bob": "left midfielder",
    "Claire": "left winger",
    "Dave": "center midfielder",
    "Eve": "striker"
}

# Swap 1: Eve and Dave trade positions
positions["Eve"], positions["Dave"] = positions["Dave"], positions["Eve"]

# Swap 2: Dave and Alice trade positions
positions["Dave"], positions["Alice"] = positions["Alice"], positions["Dave"]

# Swap 3: Claire and Bob trade positions
positions["Claire"], positions["Bob"] = positions["Bob"], positions["Claire"]

# Swap 4: Bob and Eve trade positions
positions["Bob"], positions["Eve"] = positions["Eve"], positions["Bob"]

# Swap 5: Alice and Dave trade positions
positions["Alice"], positions["Dave"] = positions["Dave"], positions["Alice"]

# Final position of Eve
print(positions["Eve"])