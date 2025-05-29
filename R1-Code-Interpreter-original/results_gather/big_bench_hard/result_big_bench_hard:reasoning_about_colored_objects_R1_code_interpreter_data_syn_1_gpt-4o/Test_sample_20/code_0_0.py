# List of objects with their colors
objects = [
    ("necklace", "black"),
    ("fidget spinner", "green"),
    ("keychain", "blue"),
    ("sheet of paper", "yellow"),
    ("stress ball", "red")
]

# Count objects that are neither yellow nor green
count = sum(1 for obj, color in objects if color not in ["yellow", "green"])

print(count)