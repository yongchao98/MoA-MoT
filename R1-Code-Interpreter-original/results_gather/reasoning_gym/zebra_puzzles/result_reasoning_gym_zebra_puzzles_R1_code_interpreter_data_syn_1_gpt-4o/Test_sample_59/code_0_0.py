# Define the houses with their attributes
houses = [
    {"smoothie": "Butterscotch", "phone": "Google Pixel 6", "name": None, "color": None},
    {"smoothie": None, "phone": None, "name": "Bob", "color": None},
    {"smoothie": "Darkness", "phone": None, "name": "Alice", "color": "Purple"},
    {"smoothie": "Dragonfruit", "phone": "OnePlus 9", "name": None, "color": "Blue"}
]

# Assign the remaining attributes based on the clues
houses[0]["name"] = "Arnold"
houses[0]["color"] = "Brown"
houses[1]["smoothie"] = "Desert"
houses[1]["phone"] = "Huawei P50"
houses[1]["color"] = "Brown"
houses[2]["phone"] = "Samsung Galaxy S21"
houses[3]["name"] = "Carol"

# Print the name of the person in the first house
print(houses[0]["name"])