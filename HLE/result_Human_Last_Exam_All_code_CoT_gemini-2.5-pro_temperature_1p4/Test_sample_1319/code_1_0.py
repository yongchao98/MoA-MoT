def find_tallest_mountain():
    # A dictionary of mountains mentioned in the Iliad and their elevations in meters.
    # Mount Olympus is included for context but will be excluded from the search.
    mountains = {
        "Mount Olympus": 2918,
        "Mount Ida (Crete)": 2456,
        "Mount Athos": 2033,
        "Mount Ida (Troy)": 1774,
        "Mount Pelion": 1624,
        "Mount Samothrace": 1611,
    }

    # Remove Mount Olympus from our consideration as per the user's request.
    mountains_without_olympus = {name: elev for name, elev in mountains.items() if name != "Mount Olympus"}

    # Find the tallest mountain among the remaining ones.
    tallest_mountain_name = max(mountains_without_olympus, key=mountains_without_olympus.get)
    tallest_mountain_elevation = mountains_without_olympus[tallest_mountain_name]

    print(f"After Mount Olympus, the tallest historical mountain mentioned in the Iliad is {tallest_mountain_name}.")
    print(f"Its elevation is {tallest_mountain_elevation} meters.")

if __name__ == "__main__":
    find_tallest_mountain()