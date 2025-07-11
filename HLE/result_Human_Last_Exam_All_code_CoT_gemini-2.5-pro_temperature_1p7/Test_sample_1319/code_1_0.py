def find_tallest_mountain():
    """
    Identifies the tallest historical mountain mentioned in the Iliad, after Mount Olympus.
    """
    # A dictionary of historical mountains from the Iliad and their elevations in meters.
    # Mount Olympus (2918m) is excluded as per the request.
    iliad_mountains = {
        'Mount Tmolus': 2159,
        'Mount Athos': 2033,
        'Mount Ida (near Troy)': 1774,
        'Mount Pelion': 1624,
        'Mount Samothrace (Fengari)': 1611
    }

    # Find the mountain with the highest elevation
    tallest_mountain_name = max(iliad_mountains, key=iliad_mountains.get)
    tallest_mountain_height = iliad_mountains[tallest_mountain_name]

    print(f"The tallest historical mountain mentioned in the Iliad after Mount Olympus is {tallest_mountain_name}.")
    print("\nHere is a comparison of their heights in meters:")

    # Sort the mountains by height in descending order for the final comparison
    sorted_mountains = sorted(iliad_mountains.items(), key=lambda item: item[1], reverse=True)

    # Create the comparison string, including each number
    comparison_string = " > ".join([f"{name} ({height}m)" for name, height in sorted_mountains])
    
    print(comparison_string)

find_tallest_mountain()