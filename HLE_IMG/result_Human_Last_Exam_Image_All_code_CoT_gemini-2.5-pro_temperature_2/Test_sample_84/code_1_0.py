def find_city_by_coat_of_arms():
    """
    Identifies a city based on the key features of its coat of arms.
    This function simulates a search through a database of heraldry.
    """

    # A small database of coats of arms with their distinct features.
    coats_of_arms_db = [
        {
            "city": "Amsterdam",
            "country": "Netherlands",
            "features": ["red shield", "black vertical stripe", "three white crosses"]
        },
        {
            "city": "Paris",
            "country": "France",
            "features": ["red shield", "ship sailing", "gold fleurs-de-lis"]
        },
        {
            "city": "Biograd na Moru",
            "country": "Croatia",
            "features": ["blue and yellow shield", "crescent moon", "crab", "olive branch", "oak branch"]
        },
        {
            "city": "Berlin",
            "country": "Germany",
            "features": ["white shield", "black bear"]
        }
    ]

    # Key features identified from the provided image.
    image_features_to_find = {"crab", "crescent moon", "blue and yellow shield"}

    # Search for a city whose coat of arms contains all the identified features.
    found_city = None
    for entry in coats_of_arms_db:
        # Check if the set of image features is a subset of the current entry's features.
        if image_features_to_find.issubset(set(entry["features"])):
            found_city = entry
            break

    if found_city:
        print(f"The features {list(image_features_to_find)} match the coat of arms of:")
        print(f"City: {found_city['city']}")
        print(f"Country: {found_city['country']}")
    else:
        print("No matching city found in the database.")

find_city_by_coat_of_arms()