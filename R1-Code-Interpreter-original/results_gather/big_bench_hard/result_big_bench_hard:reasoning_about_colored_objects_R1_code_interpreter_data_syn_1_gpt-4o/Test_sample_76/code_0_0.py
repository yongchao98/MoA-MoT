# List of items in order
items = ["yellow sunglasses", "brown fidget spinner", "silver cup", "mauve teddy bear"]

# Find the index of the silver item
silver_index = items.index("silver cup")

# Count non-yellow items to the right of the silver item
non_yellow_count = sum(1 for item in items[silver_index+1:] if "yellow" not in item)

# Print the count
print(non_yellow_count)