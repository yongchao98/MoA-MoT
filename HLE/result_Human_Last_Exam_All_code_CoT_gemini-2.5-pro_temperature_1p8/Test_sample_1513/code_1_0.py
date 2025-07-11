def find_island_town():
    """
    Finds a specific island town by filtering a list of candidates
    based on a set of geographical and naming criteria.
    """
    candidates = [
        {
            "name": "Hanga Roa",
            "island": "Easter Island",
            "distance_to_next_inhabited_km": 2075,
            "on_caldera_bay": True, # The bay is formed by the island's volcanic geography
            "bay_shares_name": True # The bay is commonly known as Hanga Roa Bay
        },
        {
            "name": "Edinburgh of the Seven Seas",
            "island": "Tristan da Cunha",
            "distance_to_next_inhabited_km": 2437,
            "on_caldera_bay": False, # Town is on a coastal plain, not a caldera bay
            "bay_shares_name": False
        },
        {
            "name": "Adamstown",
            "island": "Pitcairn Island",
            "distance_to_next_inhabited_km": 2075,
            "on_caldera_bay": False,
            "bay_shares_name": False
        },
        {
            "name": "Fira",
            "island": "Santorini",
            "distance_to_next_inhabited_km": 20, # Very close to other Greek islands
            "on_caldera_bay": True,
            "bay_shares_name": False # Town is Fira, bay is the Santorini Caldera
        }
    ]

    print("Searching for an island town based on the following criteria:")
    print("1. More than 500 kilometers from another inhabited island.")
    print("2. Sits on a bay formed by a volcanic caldera.")
    print("3. Shares the town's name with the bay.")
    print("-" * 30)

    for town in candidates:
        # Check if the town meets all the criteria
        is_remote = town["distance_to_next_inhabited_km"] > 500
        is_on_caldera_bay = town["on_caldera_bay"]
        name_matches_bay = town["bay_shares_name"]

        if is_remote and is_on_caldera_bay and name_matches_bay:
            print(f"Found a match: {town['name']}")
            print(f"Details:")
            print(f" - Distance to nearest inhabited land: {town['distance_to_next_inhabited_km']} km")
            print(f" - Sits on a bay that is part of a volcanic caldera: {town['on_caldera_bay']}")
            print(f" - The town's name is the same as the bay's name: {town['bay_shares_name']}")
            print("-" * 30)
            print(f"The final answer is: {town['name']}")
            return town['name']

    print("No town matching all criteria was found in the list.")
    return None

find_island_town()