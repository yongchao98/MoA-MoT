# Initial presents
presents = {
    "Alice": "Blue",
    "Bob": "Yellow",
    "Claire": "Black",
    "Dave": "Purple",
    "Eve": "Red"
}

# Swap 1: Bob and Claire
presents["Bob"], presents["Claire"] = presents["Claire"], presents["Bob"]

# Swap 2: Dave and Eve
presents["Dave"], presents["Eve"] = presents["Eve"], presents["Dave"]

# Swap 3: Bob and Claire
presents["Bob"], presents["Claire"] = presents["Claire"], presents["Bob"]

# Swap 4: Claire and Alice
presents["Claire"], presents["Alice"] = presents["Alice"], presents["Claire"]

# Swap 5: Eve and Dave
presents["Eve"], presents["Dave"] = presents["Dave"], presents["Eve"]

# Output the present Bob has at the end
print(presents["Bob"])