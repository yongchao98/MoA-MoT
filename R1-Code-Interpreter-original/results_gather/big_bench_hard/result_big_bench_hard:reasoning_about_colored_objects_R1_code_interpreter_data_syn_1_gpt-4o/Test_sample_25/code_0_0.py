# List of items on the floor
items = [
    ("silver", "pencil"), ("silver", "pencil"), ("silver", "pencil"),
    ("teal", "pencil"),
    ("turquoise", "mug"),
    ("yellow", "pencil"), ("yellow", "pencil"),
    ("teal", "bracelet"), ("teal", "bracelet"), ("teal", "bracelet"),
    ("silver", "mug"), ("silver", "mug"),
    ("turquoise", "jug"),
    ("turquoise", "pencil"), ("turquoise", "pencil"), ("turquoise", "pencil"),
    ("turquoise", "bracelet"), ("turquoise", "bracelet")
]

# Remove silver items
remaining_items = [item for item in items if item[0] != "silver"]

# Count the number of mugs remaining
mugs_remaining = sum(1 for item in remaining_items if item[1] == "mug")

print(mugs_remaining)