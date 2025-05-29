# Initial ball distribution
balls = {
    "Alice": "brown",
    "Bob": "purple",
    "Claire": "black",
    "Dave": "green",
    "Eve": "yellow"
}

# Swaps
# Claire and Alice swap balls
balls["Alice"], balls["Claire"] = balls["Claire"], balls["Alice"]

# Bob and Alice swap balls
balls["Alice"], balls["Bob"] = balls["Bob"], balls["Alice"]

# Eve and Dave swap balls
balls["Eve"], balls["Dave"] = balls["Dave"], balls["Eve"]

# Dave and Claire swap balls
balls["Dave"], balls["Claire"] = balls["Claire"], balls["Dave"]

# Alice and Bob swap balls
balls["Alice"], balls["Bob"] = balls["Bob"], balls["Alice"]

# Output the ball Claire has at the end
print(balls["Claire"])