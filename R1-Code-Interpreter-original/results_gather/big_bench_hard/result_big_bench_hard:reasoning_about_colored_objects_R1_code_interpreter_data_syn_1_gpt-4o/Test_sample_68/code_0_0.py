# List of items in order
items = ["gold textbook", "purple puzzle", "teal necklace", "silver pencil"]

# Find the index of the silver pencil
pencil_index = items.index("silver pencil")

# Count non-gold items to the right of the silver pencil
non_gold_items_right = [item for item in items[pencil_index+1:] if "gold" not in item]

# Output the count
print(len(non_gold_items_right))