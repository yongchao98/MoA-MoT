# Define the houses and their attributes
houses = [
    {"name": "Alice", "smoothie": "Desert", "phone": None, "flower": None},
    {"name": "Carol", "smoothie": None, "phone": "Huawei P50", "flower": "Daffodils"},
    {"name": None, "smoothie": "Darkness", "phone": "Samsung Galaxy S21", "flower": "Lilies"},
    {"name": "Arnold", "smoothie": None, "phone": None, "flower": "Carnations"}
]

# Assign the remaining attributes
houses[0]["name"] = "Bob"
houses[0]["smoothie"] = "Dragonfruit"
houses[2]["name"] = "Alice"

# Print the name of the person in House 1
print(houses[0]["name"])