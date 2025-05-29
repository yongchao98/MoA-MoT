# Define the houses and their attributes
houses = [
    {"name": "Arnold", "child": "Bella", "flower": None, "color": None},
    {"name": "Carol", "child": None, "flower": None, "color": "Brown"},
    {"name": None, "child": "Billy", "flower": None, "color": "White"},
    {"name": "Bob", "child": "Alice", "flower": "Daffodils", "color": "Purple"}
]

# Assign the remaining attributes based on the clues
houses[2]["name"] = "Alice's Mother"
houses[3]["flower"] = "Iris"
houses[3]["color"] = "Blue"

# Print the name of the person in House 1
print(houses[0]["name"])