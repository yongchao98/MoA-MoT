def find_tallest_mountain_in_iliad():
    """
    Identifies the tallest historical mountain mentioned in the Iliad, excluding Mount Olympus.
    """
    # A dictionary of mountains mentioned in the Iliad and their elevations in meters.
    mountains = {
        "Mount Olympus": 2918,
        "Mount Ida (Crete)": 2456,
        "Mount Athos": 2033,
        "Mount Ida (Troad)": 1774,
        "Mount Pelion": 1624,
        "Mount Samothrace (Fengari)": 1611
    }

    print("Comparing the elevations of historical mountains mentioned in the Iliad:")
    # Print each mountain and its height for context, fulfilling the "show the numbers" requirement.
    for name, height in mountains.items():
        print(f"- {name}: {height}m")

    # Create a new dictionary excluding Mount Olympus as per the problem's constraint.
    contenders = {name: height for name, height in mountains.items() if name != "Mount Olympus"}

    # Find the name of the mountain with the maximum height in the new dictionary.
    if not contenders:
        print("\nNo other mountains to compare.")
        return

    tallest_mountain_name = max(contenders, key=contenders.get)
    tallest_mountain_height = contenders[tallest_mountain_name]

    print(f"\nAfter excluding Mount Olympus ({mountains['Mount Olympus']}m), the next tallest is:")
    # The final equation is a comparison. We show the result of that comparison.
    print(f"Result: {tallest_mountain_name} with an elevation of {tallest_mountain_height}m.")

if __name__ == "__main__":
    find_tallest_mountain_in_iliad()