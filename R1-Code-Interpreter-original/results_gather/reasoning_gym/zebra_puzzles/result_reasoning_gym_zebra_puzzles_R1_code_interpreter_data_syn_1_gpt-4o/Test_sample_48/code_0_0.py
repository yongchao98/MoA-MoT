# Define the houses and their attributes
houses = [
    {"name": None, "color": None, "lunch": None, "nationality": None},
    {"name": "alice", "color": None, "lunch": "stir fry", "nationality": None},
    {"name": "arnold", "color": "white", "lunch": "pizza", "nationality": None},
    {"name": "bob", "color": "brown", "lunch": "soup", "nationality": "chinese"}
]

# Assign known attributes
houses[0]["nationality"] = "brit"
houses[0]["name"] = "carol"
houses[0]["color"] = "blue"
houses[0]["lunch"] = "grilled cheese"
houses[0]["nationality"] = "german"

# Print the name of the person in the first house
print(houses[0]["name"])