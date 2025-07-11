def find_tallest_iliad_mountain():
    """
    Identifies the tallest historical mountain mentioned in the Iliad, after Mount Olympus.
    """
    # A dictionary of mountains mentioned in the Iliad and their elevations in meters.
    # Mount Olympus (~2918m) is excluded as per the user's request.
    mountains = {
        "Mount Ida (Kaz Dağı)": 1774,
        "Mount Samothrace (Fengari)": 1611,
        "Mount Athos": 2033,
        "Mount Pelion": 1624,
        "Mount Ossa (Kissavos)": 1978,
        "Mount Tmolus (Bozdağ)": 2157
    }

    # Find the mountain with the highest elevation
    tallest_mountain_name = max(mountains, key=mountains.get)
    tallest_mountain_height = mountains[tallest_mountain_name]

    print(f"The mountains mentioned in the Iliad (other than Olympus) and their elevations are:")
    for name, height in mountains.items():
        print(f"- {name}: {height} meters")

    print("\nAfter comparing their heights, the tallest among them is:")
    print(f"{tallest_mountain_name} with an elevation of {tallest_mountain_height} meters.")

find_tallest_iliad_mountain()