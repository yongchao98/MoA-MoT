# Define the houses and their attributes
houses = [
    {"name": "alice", "flower": None, "lunch": "pizza", "pet": "cat"},
    {"name": "bob", "flower": None, "lunch": "grilled cheese", "pet": "fish"},
    {"name": "carol", "flower": "iris", "lunch": None, "pet": "bird"},
    {"name": "arnold", "flower": "carnations", "lunch": "soup", "pet": None}
]

# Check the conditions
conditions_met = (
    houses[0]["lunch"] == "pizza" and houses[0]["pet"] == "cat" and
    houses[1]["name"] == "bob" and
    houses[2]["name"] == "carol" and houses[2]["flower"] == "iris" and
    houses[3]["name"] == "arnold" and houses[3]["lunch"] == "soup" and houses[3]["flower"] == "carnations"
)

# Print the name of the person in house 1 if all conditions are met
if conditions_met:
    print(houses[0]["name"])
else:
    print("Conditions not met")