# Initial gifts
gifts = {
    "Alice": "blue present",
    "Bob": "green present",
    "Claire": "pink ball",
    "Dave": "yellow present",
    "Eve": "black ball"
}

# Swap 1: Bob and Dave
gifts["Bob"], gifts["Dave"] = gifts["Dave"], gifts["Bob"]

# Swap 2: Bob and Alice
gifts["Bob"], gifts["Alice"] = gifts["Alice"], gifts["Bob"]

# Swap 3: Eve and Dave
gifts["Eve"], gifts["Dave"] = gifts["Dave"], gifts["Eve"]

# Swap 4: Bob and Eve
gifts["Bob"], gifts["Eve"] = gifts["Eve"], gifts["Bob"]

# Swap 5: Eve and Claire
gifts["Eve"], gifts["Claire"] = gifts["Claire"], gifts["Eve"]

# Output the final gift Eve has
print(gifts["Eve"])