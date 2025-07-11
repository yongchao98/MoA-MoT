def find_tallest_mountain_in_iliad():
    """
    Identifies the tallest historical mountain mentioned in the Iliad, after Mount Olympus.
    """
    # A dictionary of historical mountains mentioned in the Iliad and their elevations in meters.
    # Sources for these elevations are standard encyclopedic and geographic references.
    mountains = {
        "Olympus": 2918,
        "Tmolus": 2159,
        "Athos": 2033,
        "Ida (Troad)": 1774,
        "Pelion": 1624,
        "Samothrace (Fengari)": 1611,
        "Sipylus": 1513,
    }

    # Exclude Mount Olympus as per the user's request
    other_mountains = {name: height for name, height in mountains.items() if name != "Olympus"}

    # Find the tallest mountain among the rest
    if not other_mountains:
        print("No other mountains to compare.")
        return

    tallest_mountain_name = max(other_mountains, key=other_mountains.get)
    tallest_mountain_height = other_mountains[tallest_mountain_name]

    print("After Mount Olympus, the historical mountains mentioned in the Iliad are compared below:")
    for name, height in sorted(other_mountains.items(), key=lambda item: item[1], reverse=True):
        print(f"- {name}: {height} meters")

    print("\nThe tallest among them is Mount Tmolus.")
    print("\nFinal Answer:")
    # Fulfilling the requirement to output the numbers in the final equation/result.
    print(f"Mount Tmolus = {tallest_mountain_height}")

find_tallest_mountain_in_iliad()