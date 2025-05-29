# Initial ball distribution
balls = {
    "Alice": "black",
    "Bob": "blue",
    "Claire": "orange",
    "Dave": "pink",
    "Eve": "white"
}

# Perform the swaps
# Eve and Bob swap
balls["Eve"], balls["Bob"] = balls["Bob"], balls["Eve"]

# Dave and Bob swap
balls["Dave"], balls["Bob"] = balls["Bob"], balls["Dave"]

# Dave and Alice swap
balls["Dave"], balls["Alice"] = balls["Alice"], balls["Dave"]

# Eve and Claire swap
balls["Eve"], balls["Claire"] = balls["Claire"], balls["Eve"]

# Claire and Alice swap
balls["Claire"], balls["Alice"] = balls["Alice"], balls["Claire"]

# Output the ball Eve has at the end
print(balls["Eve"])