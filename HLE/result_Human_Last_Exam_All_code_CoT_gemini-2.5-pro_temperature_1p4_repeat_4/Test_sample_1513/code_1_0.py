def solve_riddle():
    """
    This function solves the geographical riddle by analyzing the given clues.
    """
    # Clue 1: The minimum distance from another inhabited island in kilometers.
    distance_clue = 500

    # Clue 2: The town is on a bay formed by a volcanic caldera.
    # Clue 3: The town and the bay have the same name.

    # Deduction based on the clues:
    # The combination of extreme isolation (far more than 500 km), a volcanic origin,
    # and a shared name for the town and its bay points directly to one location on Earth.
    # The island is Easter Island (Rapa Nui), and its main town is Hanga Roa,
    # which is situated on Hanga Roa Bay.

    town_name = "Hanga Roa"
    island_name = "Easter Island (Rapa Nui)"

    # The final equation of the riddle is the combination of clues leading to the answer.
    # Here, we output the number from the riddle and the final answer.
    print(f"The clue specified a distance greater than {distance_clue} kilometers.")
    print(f"The island town that matches all the criteria is: {town_name}")

solve_riddle()