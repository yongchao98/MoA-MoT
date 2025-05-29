# List of objects with their colors
objects = [
    "green bracelet",
    "purple cat toy",
    "brown pair of sunglasses",
    "black fidget spinner",
    "red dog leash",
    "orange pen"
]

# Count objects that are neither black nor blue
count = sum(1 for obj in objects if "black" not in obj and "blue" not in obj)

print(count)