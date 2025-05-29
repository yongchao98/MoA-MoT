# Initial distribution
presents = {
    "Alice": "Green",
    "Bob": "Purple",
    "Claire": "Blue",
    "Dave": "Black",
    "Eve": "White"
}

# Swap 1: Eve and Bob
presents["Eve"], presents["Bob"] = presents["Bob"], presents["Eve"]

# Swap 2: Claire and Alice
presents["Claire"], presents["Alice"] = presents["Alice"], presents["Claire"]

# Swap 3: Bob and Eve
presents["Bob"], presents["Eve"] = presents["Eve"], presents["Bob"]

# Swap 4: Dave and Claire
presents["Dave"], presents["Claire"] = presents["Claire"], presents["Dave"]

# Swap 5: Alice and Eve
presents["Alice"], presents["Eve"] = presents["Eve"], presents["Alice"]

# Output the present Dave has at the end
print(presents["Dave"])