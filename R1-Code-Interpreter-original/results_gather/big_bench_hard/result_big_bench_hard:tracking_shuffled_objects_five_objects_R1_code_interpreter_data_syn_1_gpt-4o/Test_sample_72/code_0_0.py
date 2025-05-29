# Initial ball distribution
balls = {
    "Alice": "Black",
    "Bob": "Brown",
    "Claire": "Blue",
    "Dave": "Green",
    "Eve": "Yellow"
}

# Swap operations
# Eve and Claire swap
balls["Eve"], balls["Claire"] = balls["Claire"], balls["Eve"]

# Bob and Dave swap
balls["Bob"], balls["Dave"] = balls["Dave"], balls["Bob"]

# Eve and Alice swap
balls["Eve"], balls["Alice"] = balls["Alice"], balls["Eve"]

# Claire and Bob swap
balls["Claire"], balls["Bob"] = balls["Bob"], balls["Claire"]

# Alice and Eve swap
balls["Alice"], balls["Eve"] = balls["Eve"], balls["Alice"]

# Output the ball Bob has at the end
print(balls["Bob"])