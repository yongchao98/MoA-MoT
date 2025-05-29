# Define the houses and their attributes
houses = [
    {"name": "Carol", "child": "Bella", "flower": None, "color": "brown"},
    {"name": "Alice", "child": "Billy", "flower": "lilies", "color": None},
    {"name": "Arnold", "child": None, "flower": None, "color": "blue"},
    {"name": "Bob", "child": "Alice", "flower": "carnations", "color": "white"}
]

# Check the arrangement
for i, house in enumerate(houses, start=1):
    if house["name"] == "Carol":
        print(f"House {i}: {house['name']}")