# A dictionary of historical mountains mentioned in the Iliad (excluding Olympus) and their elevations in meters.
mountains = {
    "Mount Ida": 1774,
    "Mount Samothrace": 1611,
    "Mount Pelion": 1624,
    "Mount Ossa": 1978,
    "Mount Athos": 2033
}

# Find the name of the tallest mountain by finding the key with the maximum value in the dictionary.
tallest_mountain_name = max(mountains, key=mountains.get)
tallest_mountain_height = mountains[tallest_mountain_name]

print("Comparing the heights of historical mountains mentioned in the Iliad (after Olympus):")
# Loop through the dictionary and print each mountain and its height.
for name, height in mountains.items():
    print(f"- {name}: {height} meters")

print("\nThe tallest among these is:")
print(f"{tallest_mountain_name} with an elevation of {tallest_mountain_height} meters.")