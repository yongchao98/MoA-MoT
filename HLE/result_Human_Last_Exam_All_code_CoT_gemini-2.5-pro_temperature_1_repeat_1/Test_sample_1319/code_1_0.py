def find_tallest_mountain_in_iliad():
    """
    Identifies the tallest historical mountain mentioned in the Iliad, after Mount Olympus.
    """
    # Step 1 & 2: A dictionary of mountains from the Iliad (excluding Olympus) and their elevations in meters.
    mountains = {
        "Mount Ida": 1774,
        "Mount Samothrace": 1611,
        "Mount Athos": 2033,
        "Mount Pelion": 1624,
        "Mount Ossa": 1978,
    }

    # Step 3: Find the tallest mountain
    tallest_mountain_name = ""
    tallest_mountain_height = 0

    print("Comparing the elevations of mountains mentioned in the Iliad:")
    # Using max with a lambda function to find the key with the highest value
    tallest_mountain_name = max(mountains, key=mountains.get)
    tallest_mountain_height = mountains[tallest_mountain_name]

    # Step 4: Print the comparison and the result
    for name, height in mountains.items():
        print(f"- {name}: {height} meters")

    print(f"\nAfter Mount Olympus, the tallest historical mountain mentioned in the Iliad is {tallest_mountain_name} with an elevation of {tallest_mountain_height} meters.")

find_tallest_mountain_in_iliad()