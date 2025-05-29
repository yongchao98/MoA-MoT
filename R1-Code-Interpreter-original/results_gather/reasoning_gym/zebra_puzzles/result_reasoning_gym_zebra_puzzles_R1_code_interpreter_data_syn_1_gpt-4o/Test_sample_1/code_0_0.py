# Define the houses and their attributes
houses = [
    {"house": 1, "color": "white", "name": None, "cigar": None, "lunch": None},
    {"house": 2, "color": None, "name": "alice", "cigar": None, "lunch": "stir fry"},
    {"house": 3, "color": "brown", "name": "bob", "cigar": None, "lunch": None},
    {"house": 4, "color": None, "name": "arnold", "cigar": "blue master", "lunch": None}
]

# Assign the remaining attributes based on the clues
for house in houses:
    if house["house"] == 1:
        house["name"] = "carol"
    if house["house"] == 3:
        house["cigar"] = "prince"
        house["lunch"] = "grilled cheese"
    if house["house"] == 4:
        house["color"] = "blue"
        house["lunch"] = "soup"
        house["cigar"] = "dunhill"

# Print the name of the person in House 1
print(houses[0]["name"])