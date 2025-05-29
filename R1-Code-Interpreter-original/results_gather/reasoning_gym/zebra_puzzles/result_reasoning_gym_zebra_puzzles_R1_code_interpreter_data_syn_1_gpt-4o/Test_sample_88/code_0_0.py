# Define the houses and their attributes
houses = [
    {"name": None, "animal": None, "nationality": None, "color": None},
    {"name": None, "animal": None, "nationality": None, "color": None},
    {"name": None, "animal": None, "nationality": None, "color": None},
    {"name": None, "animal": None, "nationality": None, "color": None}
]

# Apply the clues
houses[2]["nationality"] = "chinese"
houses[2]["animal"] = "cat"
houses[3]["name"] = "carol"
houses[3]["animal"] = "bird"
houses[3]["color"] = "blue"
houses[1]["nationality"] = "brit"
houses[1]["color"] = "purple"
houses[0]["name"] = "arnold"
houses[0]["animal"] = "dog"
houses[1]["name"] = "alice"
houses[2]["nationality"] = "dane"
houses[2]["nationality"] = "german"

# Determine the name of the person in house 1
name_house_1 = houses[0]["name"]

print(name_house_1)