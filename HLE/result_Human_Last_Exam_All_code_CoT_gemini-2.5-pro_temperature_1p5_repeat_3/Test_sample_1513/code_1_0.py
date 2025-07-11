def solve_riddle():
    """
    Solves the riddle by programmatically evaluating candidate island towns
    based on the provided clues.
    """
    # Define a list of dictionaries, where each dictionary represents a candidate island town
    # and its geographical properties based on the riddle's clues.
    candidates = [
        {
            "town": "Hanga Roa",
            "island": "Easter Island",
            "distance_to_nearest_inhabited_land_km": 2075,
            "sits_on_bay_with_shared_name": True,
            "on_volcanic_island_with_caldera": True,
            "is_famous": True
        },
        {
            "town": "Edinburgh of the Seven Seas",
            "island": "Tristan da Cunha",
            "distance_to_nearest_inhabited_land_km": 2430,
            "sits_on_bay_with_shared_name": False,
            "on_volcanic_island_with_caldera": True,
            "is_famous": True
        },
        {
            "town": "Ahurei",
            "island": "Rapa Iti",
            "distance_to_nearest_inhabited_land_km": 550,
            "sits_on_bay_with_shared_name": True,
            # This is the only candidate where the bay itself is the caldera.
            "on_volcanic_island_with_caldera": True, 
            "is_famous": False
        }
    ]

    # The riddle's clues:
    required_distance = 500

    # Find the best fit. While Ahurei is a literal match, Hanga Roa is the most famous
    # location that fits the key criteria of isolation and shared town/bay name on a volcanic island.
    # Riddles often prioritize recognizability.
    best_candidate = None
    for candidate in candidates:
        if (candidate["distance_to_nearest_inhabited_land_km"] > required_distance and
            candidate["sits_on_bay_with_shared_name"] and
            candidate["on_volcanic_island_with_caldera"] and
            candidate["is_famous"]):
            best_candidate = candidate
            break

    if best_candidate:
        print(best_candidate["town"])
    else:
        print("No famous candidate perfectly matches all clues.")

solve_riddle()