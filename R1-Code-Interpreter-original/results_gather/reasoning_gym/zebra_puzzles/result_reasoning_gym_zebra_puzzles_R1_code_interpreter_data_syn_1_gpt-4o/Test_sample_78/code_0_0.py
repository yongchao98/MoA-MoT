# Define the houses and their attributes
houses = [
    {"name": "arnold", "nationality": "brit", "phone": None, "cigar": "dunhill"},
    {"name": None, "nationality": "dane", "phone": "samsung galaxy s21", "cigar": None},
    {"name": None, "nationality": "chinese", "phone": "huawei p50", "cigar": "blue master"},
    {"name": "bob", "nationality": None, "phone": "oneplus 9", "cigar": None}
]

# Assign Alice to the house where Bob is
houses[3]["name"] = "alice"

# Assign the remaining attributes based on the clues
houses[0]["phone"] = None
houses[0]["cigar"] = "dunhill"
houses[1]["cigar"] = None
houses[2]["cigar"] = "blue master"
houses[3]["cigar"] = None

# Print the name of the person in House 1
print(houses[0]["name"])