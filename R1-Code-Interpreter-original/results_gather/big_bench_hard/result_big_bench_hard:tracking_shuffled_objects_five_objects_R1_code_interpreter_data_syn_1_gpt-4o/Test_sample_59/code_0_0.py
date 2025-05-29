# Initial presents
presents = {
    "Alice": "Brown",
    "Bob": "Black",
    "Claire": "Green",
    "Dave": "Purple",
    "Eve": "Yellow"
}

# Swap 1: Bob and Dave
presents["Bob"], presents["Dave"] = presents["Dave"], presents["Bob"]

# Swap 2: Bob and Eve
presents["Bob"], presents["Eve"] = presents["Eve"], presents["Bob"]

# Swap 3: Dave and Claire
presents["Dave"], presents["Claire"] = presents["Claire"], presents["Dave"]

# Swap 4: Dave and Bob
presents["Dave"], presents["Bob"] = presents["Bob"], presents["Dave"]

# Swap 5: Alice and Dave
presents["Alice"], presents["Dave"] = presents["Dave"], presents["Alice"]

# Output the present Dave has at the end
print(presents["Dave"])