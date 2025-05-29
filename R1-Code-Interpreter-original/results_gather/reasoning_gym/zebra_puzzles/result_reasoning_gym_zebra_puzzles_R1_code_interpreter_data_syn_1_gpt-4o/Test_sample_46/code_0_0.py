# Define the houses and their attributes
houses = [
    {"house": 1, "name": None, "child": "Alice", "phone": "OnePlus 9", "color": None},
    {"house": 2, "name": "Carol", "child": "Timothy", "phone": "Google Pixel 6", "color": "blue"},
    {"house": 3, "name": "Bob", "child": None, "phone": None, "color": "purple"},
    {"house": 4, "name": None, "child": "Billy", "phone": "Huawei P50", "color": "brown"}
]

# Assign the remaining names
for house in houses:
    if house["name"] is None:
        if house["house"] == 1:
            house["name"] = "Arnold"
        elif house["house"] == 4:
            house["name"] = "Alice"

# Print the name of the person in House 1
print(houses[0]["name"])