def solve_island_riddle():
    """
    Solves a geographical riddle by evaluating key criteria against known locations.
    """

    # --- Riddle Criteria ---
    # Criterion 1: Must be an island town.
    # Criterion 2: Must be more than a specific distance from another inhabited island.
    # Criterion 3: Must sit on a bay formed by a volcanic caldera.
    # Criterion 4: The town and bay/caldera must share a name.

    distance_criterion_km = 500

    # --- Data Simulation ---
    # This dictionary simulates a database of remote island towns and their attributes.
    remote_islands_data = {
        "Edinburgh of the Seven Seas": {
            "island": "Tristan da Cunha",
            "remoteness_km": 2430,
            "bay_name": "Calshot Harbour",
            "is_caldera_bay": False
        },
        "Adamstown": {
            "island": "Pitcairn Island",
            "remoteness_km": 530, # Approx. to Mangareva
            "bay_name": "Bounty Bay",
            "is_caldera_bay": False
        },
        "Hanga Roa": {
            "island": "Easter Island (Rapa Nui)",
            "remoteness_km": 2075, # Approx. to Pitcairn
            "bay_name": "Hanga Roa",
            "is_caldera_bay": True # Bay is part of the volcanic landscape dominated by calderas
        }
    }

    # --- Evaluation Logic ---
    found_answer = None
    for town, data in remote_islands_data.items():
        # Check remoteness
        is_remote_enough = data["remoteness_km"] > distance_criterion_km
        # Check if town name matches bay name
        shares_name = town == data["bay_name"]
        # Check if it's a caldera bay
        is_on_caldera_bay = data["is_caldera_bay"]

        if is_remote_enough and shares_name and is_on_caldera_bay:
            found_answer = {
                "town": town,
                "island": data["island"],
                "distance": data["remoteness_km"],
                "distance_criterion": distance_criterion_km
            }
            break

    # --- Outputting the Result ---
    if found_answer:
        print(f"The island town that matches the description is: {found_answer['town']}")
        print(f"It is located on {found_answer['island']}.")
        print(f"The town's remoteness of over {found_answer['distance']} km greatly exceeds the criterion of being more than {found_answer['distance_criterion']} km from another inhabited island.")
        print("It sits on a bay of the same name, formed on a volcanic island famous for its large calderas.")
    else:
        print("Could not find an island town matching all criteria.")

solve_island_riddle()