# Initial gifts
gifts = {
    "Alice": "green present",
    "Bob": "brown present",
    "Claire": "white present",
    "Dave": "black ball",
    "Eve": "orange ball"
}

# Swaps
# Dave and Claire swap
gifts["Dave"], gifts["Claire"] = gifts["Claire"], gifts["Dave"]

# Dave and Eve swap
gifts["Dave"], gifts["Eve"] = gifts["Eve"], gifts["Dave"]

# Alice and Claire swap
gifts["Alice"], gifts["Claire"] = gifts["Claire"], gifts["Alice"]

# Eve and Alice swap
gifts["Eve"], gifts["Alice"] = gifts["Alice"], gifts["Eve"]

# Bob and Dave swap
gifts["Bob"], gifts["Dave"] = gifts["Dave"], gifts["Bob"]

# Output the gift Eve has
print(gifts["Eve"])