# Define the houses and their attributes
houses = [
    {"name": "Arnold", "color": "white", "cigar": "Dunhill", "child": None},
    {"name": "Alice", "color": "brown", "cigar": "Pall Mall", "child": "Timothy"},
    {"name": "Bob", "color": "blue", "cigar": None, "child": "Billy"},
    {"name": None, "color": "purple", "cigar": "Blue Master", "child": None}
]

# Find the name of the person in House 1
house_1_name = houses[0]["name"]
print(house_1_name)