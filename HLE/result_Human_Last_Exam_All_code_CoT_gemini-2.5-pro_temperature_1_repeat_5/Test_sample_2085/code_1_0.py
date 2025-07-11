def find_french_equivalent():
    """
    Finds the modern French equivalent of a historical name by looking it up
    in a predefined dictionary.
    """
    # This dictionary simulates a database of historical name equivalents.
    # Odobeccus is the Latinized form of a Frankish name.
    historical_name_map = {
        "Odobeccus": "Eudes"
    }

    gaulish_name = "Odobeccus"

    # Retrieve the modern French name from the map.
    french_name = historical_name_map.get(gaulish_name, "Unknown")

    if french_name != "Unknown":
        print(f"The French equivalent of the name '{gaulish_name}' is '{french_name}'.")
    else:
        print(f"Could not find a French equivalent for '{gaulish_name}'.")

find_french_equivalent()