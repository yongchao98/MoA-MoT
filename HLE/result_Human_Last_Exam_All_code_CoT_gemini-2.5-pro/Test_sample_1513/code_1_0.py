def solve_riddle():
    """
    This function solves the geographical riddle based on the provided clues.
    """
    # Clue 1: Extreme remoteness (> 500km from another inhabited island).
    # This points to Easter Island (Rapa Nui).
    # The distance from Easter Island to the nearest inhabited land (Pitcairn Islands)
    # is over 2,000 kilometers.
    distance_to_nearest_inhabited_land_km = 2075 # Approx. to Pitcairn

    # Clue 2: Sits on a bay formed by a volcanic caldera.
    # Easter Island is a volcanic island.

    # Clue 3: The town, bay, and caldera share the town's name.
    # The main town on Easter Island is Hanga Roa.
    # The town is located on Hanga Roa Bay.
    # The name matches for the town and the bay.
    town_name = "Hanga Roa"
    bay_name = "Hanga Roa Bay"

    print(f"The riddle describes an island town more than {500} kilometers from another inhabited island.")
    print(f"The town sits on a bay that shares its name.")
    print(f"Based on these clues, the answer is the main town on Easter Island.")
    print(f"Town Name: {town_name}")

solve_riddle()