def find_tallest_mountain():
    # A dictionary of historical mountains mentioned in the Iliad and their elevations in meters.
    # Mount Olympus (2918m) is excluded as per the user's request.
    # Sources for mountains mentioned:
    # - Mount Ida: Zeus watches the battle from its peak.
    # - Mount Samothrace (Fengari): Poseidon watches the battle from here.
    # - Mount Athos: Mentioned in a simile by Hera.
    mountains = {
        "Mount Ida": 1774,
        "Mount Samothrace (Fengari)": 1611,
        "Mount Athos": 2033
    }

    # Find the mountain with the highest elevation
    if not mountains:
        print("No mountains to compare.")
        return

    tallest_mountain_name = max(mountains, key=mountains.get)
    tallest_mountain_height = mountains[tallest_mountain_name]

    print("After Mount Olympus, the historical mountains mentioned in the Iliad are:")
    for name, height in mountains.items():
        print(f"- {name}: {height} meters")

    print("\nComparing their heights:")
    # The "final equation" is finding the maximum value among the heights.
    # We will show the numbers being compared.
    heights = list(mountains.values())
    print(f"max({', '.join(map(str, heights))}) = {tallest_mountain_height}")


    print(f"\nThe tallest historical mountain mentioned in the Iliad, after Mount Olympus, is {tallest_mountain_name} with an elevation of {tallest_mountain_height} meters.")

find_tallest_mountain()