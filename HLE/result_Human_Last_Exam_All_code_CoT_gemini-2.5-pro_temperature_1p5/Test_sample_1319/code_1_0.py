def find_tallest_mountain_in_iliad():
    """
    Identifies the tallest historical mountain mentioned in the Iliad, excluding Mount Olympus.
    The function stores the mountains and their elevations in a dictionary, then iterates
    through them to find the highest peak after Olympus.
    """
    # A dictionary of prominent historical mountains mentioned in the Iliad and their elevations in meters.
    mountains = {
        "Mount Olympus": 2918,
        "Mount Athos": 2033,
        "Mount Ida": 1774,
        "Mount Pelion": 1624,
        "Mount Samothrace": 1611
    }

    # Variables to store the result
    tallest_after_olympus_name = ""
    tallest_after_olympus_height = 0

    # Find the tallest mountain, skipping Olympus
    for name, height in mountains.items():
        if name != "Mount Olympus":
            if height > tallest_after_olympus_height:
                tallest_after_olympus_height = height
                tallest_after_olympus_name = name

    # Print the result and the final numbers used
    print(f"After Mount Olympus, the tallest historical mountain mentioned in the Iliad is {tallest_after_olympus_name}.")
    print(f"The comparison shows that the height of {tallest_after_olympus_name} is {tallest_after_olympus_height} meters.")
    print(f"This makes it taller than Mount Ida ({mountains['Mount Ida']}m), Mount Pelion ({mountains['Mount Pelion']}m), and Mount Samothrace ({mountains['Mount Samothrace']}m).")

if __name__ == '__main__':
    find_tallest_mountain_in_iliad()