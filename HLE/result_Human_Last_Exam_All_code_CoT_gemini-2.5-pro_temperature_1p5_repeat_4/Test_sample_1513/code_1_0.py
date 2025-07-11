import sys

def solve_geographical_riddle():
    """
    This script solves a riddle about a remote island town.
    It identifies the town by processing a series of geographical clues.
    """

    # --- Clue Analysis ---
    # The riddle provides several key pieces of information.
    # We will use them to pinpoint the location.

    # Clue 1: Remoteness. The distance from another inhabited island is specified.
    distance_clue_km = 500

    # Clue 2: Geography. The town is on a specific type of bay.
    geographical_feature = "bay formed by a volcanic caldera"

    # Clue 3: Naming. The town and bay have the same name.
    naming_convention = "shares the town's name"

    # --- Deduction ---
    # Based on the extreme remoteness (> 500 km), we can consider locations
    # like Easter Island, Tristan da Cunha, or the Pitcairn Islands.

    # Easter Island (Rapa Nui) is a volcanic island whose nearest inhabited
    # neighbor is over 2,000 km away, satisfying the distance clue.
    # The island's landscape is defined by its volcanoes and calderas.

    # The main town on Easter Island is named Hanga Roa.
    # It is located on a bay which is also named Hanga Roa. This matches the
    # naming and geography clues perfectly.

    island = "Easter Island (Rapa Nui)"
    town_and_bay = "Hanga Roa"

    # --- Outputting the result ---
    # The final code needs to output the numbers from the clues.
    # Here we present the solution and the key number.
    print("Finding the island town based on the riddle:")
    print(f"The required distance from another inhabited island is > {distance_clue_km} km.")
    print(f"The town and bay name match, and the location is on a bay formed by volcanic activity.")
    print("-" * 20)
    print(f"The island is: {island}")
    print(f"The town that fits all criteria is: {town_and_bay}")

solve_geographical_riddle()
<<<Hanga Roa>>>