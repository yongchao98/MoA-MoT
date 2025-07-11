def find_tallest_mountain_in_iliad():
    """
    Identifies the tallest historical mountain mentioned in the Iliad,
    after Mount Olympus, and prints the result.
    """
    # Data of mountains mentioned in the Iliad and their elevations in meters.
    # Source: Encyclopedia Britannica, Wikipedia, and other geographical sources.
    mountains = {
        "Mount Olympus": 2918,
        "Mount Parnassus": 2457,
        "Mount Tmolus": 2157,
        "Mount Athos": 2033,
        "Mount Ossa": 1978,
        "Mount Ida (Turkey)": 1774,
        "Mount Pelion": 1624,
        "Mount Samothrace": 1611,
    }

    # The problem asks for the tallest mountain *after* Mount Olympus.
    # So, we remove it from our dictionary before comparison.
    if "Mount Olympus" in mountains:
        del mountains["Mount Olympus"]

    # Find the mountain with the highest elevation from the remaining list.
    # The max() function with a key lambda function is perfect for finding the
    # key in a dictionary corresponding to the maximum value.
    tallest_mountain_name = max(mountains, key=mountains.get)
    tallest_mountain_height = mountains[tallest_mountain_name]

    print(f"After Mount Olympus, the tallest historical mountain mentioned in the Iliad is {tallest_mountain_name}.")
    print(f"It has an elevation of {tallest_mountain_height} meters.")

find_tallest_mountain_in_iliad()