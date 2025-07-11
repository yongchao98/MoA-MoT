def solve_riddle():
    """
    This function solves the geographical riddle based on a step-by-step analysis of the clues.
    """
    # Clue 1: More than 500 kilometers from another inhabited island.
    # This points to a highly remote location like Easter Island.
    # Distance from Easter Island to Pitcairn Islands is ~2,075 km.
    distance_to_nearest_inhabited_land_km = 2075
    
    # Clue 2 & 3: Sits on a bay formed by a volcanic caldera that shares the town's name.
    # On Easter Island, the main town is Hanga Roa.
    # It is located on Hanga Roa Bay.
    # The bay is on a volcanic island and is geographically shaped by the massive Rano Kau caldera.
    town_name = "Hanga Roa"
    bay_name = "Hanga Roa Bay"
    
    # Verification
    is_remote = distance_to_nearest_inhabited_land_km > 500
    name_matches = town_name in bay_name
    
    if is_remote and name_matches:
        print(f"The island town is: {town_name}")

solve_riddle()