# Initial ball distribution
balls = {
    "Alice": "Purple",
    "Bob": "Brown",
    "Claire": "Blue",
    "Dave": "Pink",
    "Eve": "Yellow",
    "Fred": "White",
    "Gertrude": "Black"
}

# Perform the swaps
# Alice and Eve swap
balls["Alice"], balls["Eve"] = balls["Eve"], balls["Alice"]

# Dave and Fred swap
balls["Dave"], balls["Fred"] = balls["Fred"], balls["Dave"]

# Alice and Gertrude swap
balls["Alice"], balls["Gertrude"] = balls["Gertrude"], balls["Alice"]

# Fred and Claire swap
balls["Fred"], balls["Claire"] = balls["Claire"], balls["Fred"]

# Bob and Fred swap
balls["Bob"], balls["Fred"] = balls["Fred"], balls["Bob"]

# Gertrude and Bob swap
balls["Gertrude"], balls["Bob"] = balls["Bob"], balls["Gertrude"]

# Gertrude and Claire swap
balls["Gertrude"], balls["Claire"] = balls["Claire"], balls["Gertrude"]

# Output the ball Fred has at the end
print(balls["Fred"])