def find_tallest_mountain_after_olympus():
    """
    Identifies the tallest historical mountain mentioned in the Iliad, after Mount Olympus.
    """
    # A dictionary of mountains mentioned in the Iliad and their approximate elevations in meters.
    mountains = {
        "Mount Olympus": 2918,
        "Mount Athos": 2033,
        "Mount Ida (Kaz Dağı)": 1774,
        "Mount Samothrace (Fengari)": 1611
    }

    # Remove Mount Olympus to find the next tallest
    olympus_height = mountains.pop("Mount Olympus")
    print(f"First, we consider the highest peak, Mount Olympus, at {olympus_height} meters.")
    print("Next, we compare the other mountains mentioned in the Iliad:")
    print("-" * 30)

    # Find the tallest mountain among the rest
    tallest_mountain_name = ""
    tallest_mountain_height = 0
    for name, height in mountains.items():
        print(f"{name}: {height} meters")
        if height > tallest_mountain_height:
            tallest_mountain_height = height
            tallest_mountain_name = name

    print("-" * 30)
    print(f"Comparing the heights {mountains['Mount Athos']}, {mountains['Mount Ida (Kaz Dağı)']}, and {mountains['Mount Samothrace (Fengari)']}, the highest is {tallest_mountain_height}.")
    print(f"\nTherefore, after Mount Olympus, the tallest historical mountain mentioned in the Iliad is {tallest_mountain_name}.")

find_tallest_mountain_after_olympus()