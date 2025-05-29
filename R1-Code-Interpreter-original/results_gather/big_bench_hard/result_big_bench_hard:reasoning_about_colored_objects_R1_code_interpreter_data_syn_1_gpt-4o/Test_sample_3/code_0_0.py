# List of items
items = ["fuchsia teddy bear", "red jug", "gold puzzle", "burgundy bracelet", "green notebook"]

# Find the index of the red item
red_index = items.index("red jug")

# Count non-magenta items to the right of the red item
non_magenta_count = 0
for item in items[red_index + 1:]:
    if "magenta" not in item:
        non_magenta_count += 1

# Output the count
print(non_magenta_count)