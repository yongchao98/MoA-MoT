def find_tallest_iliad_mountain():
    """
    Identifies the tallest historical mountain mentioned in the Iliad, after Mount Olympus.
    """
    # Mountains explicitly mentioned in the Iliad as locations (excluding Olympus)
    # and their respective elevations in meters.
    mountains = {
        "Mount Athos": 2033,
        "Mount Ida (near Troy)": 1774,
        "Mount Pelion": 1624,
        "Mount Samothrace": 1611
    }

    # Find the tallest mountain from the dictionary
    tallest_mountain_name = ""
    max_elevation = 0
    for name, elevation in mountains.items():
        if elevation > max_elevation:
            max_elevation = elevation
            tallest_mountain_name = name

    # Print the explanation and the result
    print("After Mount Olympus (2918 m), the historical mountains mentioned by name in the Iliad are:")
    for name, elevation in mountains.items():
        print(f"- {name}: {elevation} m")
    
    print("\nBy comparing their elevations, we can determine the tallest.")
    print(f"The tallest among these is {tallest_mountain_name} with an elevation of {max_elevation} meters.")

if __name__ == "__main__":
    find_tallest_iliad_mountain()
