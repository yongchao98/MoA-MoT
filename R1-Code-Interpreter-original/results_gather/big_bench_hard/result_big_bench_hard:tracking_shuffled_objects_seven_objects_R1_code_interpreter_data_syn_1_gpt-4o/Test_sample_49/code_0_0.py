# Initial gifts
gifts = {
    "Alice": "white present",
    "Bob": "orange ball",
    "Claire": "black ball",
    "Dave": "yellow present",
    "Eve": "green present",
    "Fred": "purple present",
    "Gertrude": "pink ball"
}

# Swaps
gifts["Fred"], gifts["Alice"] = gifts["Alice"], gifts["Fred"]
gifts["Eve"], gifts["Claire"] = gifts["Claire"], gifts["Eve"]
gifts["Dave"], gifts["Bob"] = gifts["Bob"], gifts["Dave"]
gifts["Gertrude"], gifts["Alice"] = gifts["Alice"], gifts["Gertrude"]
gifts["Eve"], gifts["Bob"] = gifts["Bob"], gifts["Eve"]
gifts["Eve"], gifts["Gertrude"] = gifts["Gertrude"], gifts["Eve"]
gifts["Dave"], gifts["Bob"] = gifts["Bob"], gifts["Dave"]

# Final gift with Bob
print(gifts["Bob"])