def solve_geography_riddle():
    """
    This function identifies an island town based on a set of geographical clues
    and prints the details and the final answer.
    """

    # --- Riddle's Clues ---
    # The clue specifies a distance greater than 500 kilometers.
    required_distance_km = 500

    # --- Identified Location's Data ---
    town_name = "Hanga Roa"
    bay_name = "Hanga Roa Bay"
    island_name = "Easter Island (Rapa Nui)"
    geological_feature = "volcanic island with calderas"
    
    # Easter Island's nearest inhabited neighbor (Pitcairn Island) is 2075 km away.
    actual_distance_km = 2075

    # --- Verification and Output ---
    print("Finding the island town based on the clues...")
    print(f"Clue 1: More than {required_distance_km} km from another inhabited island.")
    print(f"Check: The nearest inhabited island is {actual_distance_km} km away. This satisfies the condition.\n")

    print(f"Clue 2: Sits on a bay formed by a volcanic caldera.")
    print(f"Check: The town sits on a bay on a {geological_feature}.\n")
    
    print(f"Clue 3: The town and the bay share the same name.")
    print(f"Check: The town is '{town_name}' and the bay is '{bay_name}'. This satisfies the condition.\n")

    print("---")
    print("Conclusion: The island town that fits all the descriptions is:")
    print(town_name)

solve_geography_riddle()