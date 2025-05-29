# Define the houses and their attributes
houses = [
    {"name": None, "lunch": None, "smoothie": None, "cigar": None},
    {"name": None, "lunch": None, "smoothie": None, "cigar": None},
    {"name": None, "lunch": None, "smoothie": None, "cigar": None},
    {"name": None, "lunch": None, "smoothie": None, "cigar": None}
]

# Apply the clues
houses[3]["lunch"] = "pizza"
houses[3]["smoothie"] = "darkness"
houses[2]["name"] = "arnold"
houses[1]["cigar"] = "dunhill"
houses[1]["name"] = "alice"
houses[0]["cigar"] = "blue master"
houses[0]["name"] = "bob"
houses[0]["lunch"] = "soup"
houses[0]["smoothie"] = "butterscotch"
houses[2]["cigar"] = "prince"
houses[2]["lunch"] = "grilled cheese"

# Print the name of the person in House 1
print(houses[0]["name"])