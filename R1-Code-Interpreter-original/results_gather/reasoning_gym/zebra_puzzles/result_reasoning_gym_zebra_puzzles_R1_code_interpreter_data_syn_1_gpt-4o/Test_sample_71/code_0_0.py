# Define the houses and their attributes
houses = [
    {"name": None, "smoothie": None, "flower": None, "animal": None},
    {"name": None, "smoothie": None, "flower": None, "animal": None},
    {"name": None, "smoothie": None, "flower": None, "animal": None},
    {"name": None, "smoothie": None, "flower": None, "animal": None}
]

# Apply the clues
houses[2]["smoothie"] = "Dragonfruit"
houses[2]["animal"] = "cat"
houses[3]["smoothie"] = "Desert"
houses[3]["animal"] = "dog"
houses[1]["name"] = "Bob"
houses[1]["smoothie"] = "Butterscotch"
houses[1]["flower"] = "carnations"
houses[1]["animal"] = "bird"
houses[0]["name"] = "Carol"
houses[0]["smoothie"] = "Darkness"
houses[0]["name"] = "Alice"
houses[0]["flower"] = "iris"
houses[2]["flower"] = "lilies"

# Check the name of the person in house 1
print(houses[0]["name"])