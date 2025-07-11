def find_island_town():
    """
    Finds a remote island town based on specific geographical criteria by filtering a list of candidates.
    """
    # A list of dictionaries representing candidate island towns and their properties.
    # The 'bay_is_caldera' key is a simplification representing if the bay is formed by or directly associated with a caldera.
    candidates = [
        {
            'town_name': 'Hanga Roa',
            'remoteness_km': 2075,
            'bay_shares_name': True,
            'bay_is_caldera': True,
            'details': 'Town is Hanga Roa, bay is Hanga Roa Bay. The island is volcanic and the town lies near the Rano Kau caldera.'
        },
        {
            'town_name': 'Edinburgh of the Seven Seas',
            'remoteness_km': 2430,
            'bay_shares_name': False,
            'bay_is_caldera': False,
            'details': 'Town sits on a volcanic plain on Tristan da Cunha, not in a bay sharing its name.'
        },
        {
            'town_name': 'Adamstown',
            'remoteness_km': 2075, # Distance to Easter Island
            'bay_shares_name': False,
            'bay_is_caldera': False,
            'details': 'Located on Pitcairn Island, the main landing is Bounty Bay.'
        },
        {
            'town_name': 'Fira',
            'remoteness_km': 130, # Distance to Crete
            'bay_shares_name': False,
            'bay_is_caldera': True,
            'details': 'Sits on the edge of the Santorini caldera, but is not over 500km from another inhabited island.'
        },
         {
            'town_name': 'Jamestown',
            'remoteness_km': 1900, # Distance to Africa
            'bay_shares_name': False,
            'bay_is_caldera': False,
            'details': 'Sits in James Bay on St. Helena, which is not a caldera bay.'
        }
    ]

    print("Searching for the island town based on the following criteria:")
    print("1. Remoteness > 500 km")
    print("2. Bay is formed by a volcanic caldera")
    print("3. Bay shares the town's name")
    print("-" * 30)

    found_town = None
    for town in candidates:
        # Check if all conditions are met
        is_remote = town['remoteness_km'] > 500
        is_caldera_bay = town['bay_is_caldera']
        is_name_shared = town['bay_shares_name']

        if is_remote and is_caldera_bay and is_name_shared:
            found_town = town
            break
            
    if found_town:
        print(f"Match found: {found_town['town_name']}")
        print(f"Details: {found_town['details']}")
    else:
        print("No town in the list matches all criteria.")

find_island_town()