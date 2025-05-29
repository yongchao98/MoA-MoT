# Define the constraints
constraints = {
    "mangoes": 3,
    "plums": [4, 5, 6, 7],
    "watermelons": [5, 6, 7],
    "oranges": [6, 7],
    "kiwis": [7],
    "apples": [1, 2, 4, 5, 6],
    "loquats": [7]
}

# Initialize the order list
order = [None] * 7

# Place mangoes in the third position
order[2] = "mangoes"

# Place loquats in the seventh position
order[6] = "loquats"

# Place kiwis in the sixth position
order[5] = "kiwis"

# Place oranges in the fifth position
order[4] = "oranges"

# Place watermelons in the fourth position
order[3] = "watermelons"

# Place plums in the second position
order[1] = "plums"

# Place apples in the first position
order[0] = "apples"

# Check which fruit is the second-cheapest
second_cheapest = order[1]

print(second_cheapest)