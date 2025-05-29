# Initial gifts
gifts = {
    "Alice": "Brown",
    "Bob": "Black",
    "Claire": "Green",
    "Dave": "Purple",
    "Eve": "Yellow"
}

# Swaps
gifts["Bob"], gifts["Dave"] = gifts["Dave"], gifts["Bob"]
gifts["Bob"], gifts["Eve"] = gifts["Eve"], gifts["Bob"]
gifts["Dave"], gifts["Claire"] = gifts["Claire"], gifts["Dave"]
gifts["Dave"], gifts["Bob"] = gifts["Bob"], gifts["Dave"]
gifts["Alice"], gifts["Dave"] = gifts["Dave"], gifts["Alice"]

# Output the final gift Bob has
print(gifts["Bob"])