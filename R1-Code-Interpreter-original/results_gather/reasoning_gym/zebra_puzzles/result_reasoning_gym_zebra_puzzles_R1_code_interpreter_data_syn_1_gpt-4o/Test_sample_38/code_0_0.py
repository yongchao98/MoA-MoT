# Initialize the houses with empty values
houses = [
    {"name": None, "lunch": None, "child": None, "flower": None},
    {"name": None, "lunch": None, "child": None, "flower": None},
    {"name": None, "lunch": None, "child": None, "flower": None},
    {"name": None, "lunch": None, "child": None, "flower": None}
]

# Apply the clues
houses[0]["flower"] = "iris"  # Clue 9
houses[0]["name"] = "arnold"  # Clue 2, 11, and 3
houses[0]["lunch"] = "pizza"  # Clue 3
houses[0]["child"] = "billy"  # Clue 2 and 11

houses[2]["child"] = "alice"  # Clue 1

houses[3]["lunch"] = "soup"  # Clue 8
houses[3]["child"] = "timothy"  # Clue 4

# Determine the remaining attributes
for i, house in enumerate(houses):
    if house["name"] is None:
        if "carol" not in [h["name"] for h in houses]:
            house["name"] = "carol"
            house["flower"] = "daffodils"  # Clue 10
        else:
            house["name"] = "bob"
    if house["lunch"] is None:
        if "grilled cheese" not in [h["lunch"] for h in houses]:
            house["lunch"] = "grilled cheese"  # Clue 6
        else:
            house["lunch"] = "stir fry"  # Clue 5
            house["child"] = "bella"  # Clue 5

# Print the name of the person in the first house
print(houses[0]["name"])