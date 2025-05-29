# Define the houses and their attributes
houses = [
    {"name": None, "flower": None, "lunch": None, "cigar": None},
    {"name": None, "flower": None, "lunch": "stir fry", "cigar": "Pall Mall"},
    {"name": "Carol", "flower": "daffodils", "lunch": None, "cigar": "Blue Master"},
    {"name": None, "flower": "lilies", "lunch": "soup", "cigar": None}
]

# Assign known values based on clues
houses[0]["name"] = "Arnold"
houses[0]["lunch"] = "pizza"

# Check the arrangement
print(houses[0]["name"])