def find_tallest_iliad_mountain():
    """
    Identifies the tallest historical mountain mentioned in the Iliad, excluding Mount Olympus.
    """
    # Mountains mentioned in the Iliad with their real-world elevations in meters.
    # Mount Olympus is included for context but will be excluded from the main comparison.
    mountains = {
        "Mount Olympus": 2918,
        "Mount Athos": 2033,
        "Mount Ida (Turkey)": 1774,
        "Mount Othrys": 1726,
        "Mount Pelion": 1624,
        "Samothrace (Fengari Peak)": 1611,
    }

    # The problem asks to find the tallest mountain *after* Mount Olympus.
    # We create a new dictionary without Olympus to perform the comparison.
    contender_mountains = {name: height for name, height in mountains.items() if name != "Mount Olympus"}

    print("Comparing the heights of historical mountains mentioned in the Iliad (excluding Mount Olympus):")
    for name, height in contender_mountains.items():
        print(f"- {name}: {height} meters")

    # Find the mountain with the maximum height from the contenders.
    tallest_mountain_name = max(contender_mountains, key=contender_mountains.get)
    tallest_mountain_height = contender_mountains[tallest_mountain_name]

    print("\n---")
    print("The tallest mountain among the contenders is:")
    print(f"Mountain: {tallest_mountain_name}")
    print(f"Elevation: {tallest_mountain_height} meters")
    print(f"\nFinal Equation: Max({contender_mountains['Mount Athos']}, {contender_mountains['Mount Ida (Turkey)']}, {contender_mountains['Mount Othrys']}, {contender_mountains['Mount Pelion']}, {contender_mountains['Samothrace (Fengari Peak)']}) = {tallest_mountain_height}")

find_tallest_iliad_mountain()