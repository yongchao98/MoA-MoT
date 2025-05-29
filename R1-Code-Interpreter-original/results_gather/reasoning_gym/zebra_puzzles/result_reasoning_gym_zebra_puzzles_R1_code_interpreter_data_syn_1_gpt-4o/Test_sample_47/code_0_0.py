# Define the houses and their attributes
houses = [
    {"name": None, "color": None, "drink": None, "smoothie": "Darkness"},
    {"name": "Bob", "color": "Purple", "drink": "Water", "smoothie": "Dragonfruit"},
    {"name": None, "color": "White", "drink": "Milk", "smoothie": "Desert"},
    {"name": "Arnold", "color": "Blue", "drink": "Boba Tea", "smoothie": None}
]

# Assign the remaining attributes based on the clues
houses[0]["name"] = "Carol"

# Print the name of the person in House 1
print(houses[0]["name"])