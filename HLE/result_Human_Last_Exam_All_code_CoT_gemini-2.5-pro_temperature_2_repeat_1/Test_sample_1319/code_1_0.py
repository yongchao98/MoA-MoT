# A dictionary mapping the historical mountains mentioned in the Iliad (besides Olympus)
# to their respective peak elevations in meters.
mountains = {
    "Mount Ida (Troad)": 1774,
    "Mount Samothrace (Fengari)": 1611,
    "Mount Athos": 2033,
    "Mount Pelion": 1624
}

# Find the name of the tallest mountain by finding the key with the maximum value.
tallest_mountain_name = max(mountains, key=mountains.get)
tallest_mountain_elevation = mountains[tallest_mountain_name]

print("Comparing the elevations of mountains mentioned in the Iliad (after Olympus):")
# Print each mountain and its elevation.
for name, elevation in mountains.items():
    print(f"- {name}: {elevation} meters")

print("\nThe comparison shows the following order from tallest to shortest:")
# Sort the mountains by elevation in descending order for a clear comparison.
sorted_mountains = sorted(mountains.items(), key=lambda item: item[1], reverse=True)
for i, (name, elevation) in enumerate(sorted_mountains):
    # This loop outputs each number in the "final equation" or comparison.
    print(f"{i+1}. {name} at {elevation}m")


print(f"\nAfter Mount Olympus, the tallest historical mountain mentioned in the Iliad is {tallest_mountain_name} with an elevation of {tallest_mountain_elevation} meters.")