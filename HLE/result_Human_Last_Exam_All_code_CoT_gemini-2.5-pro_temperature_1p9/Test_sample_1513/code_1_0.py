def solve_geography_riddle():
    """
    This function solves the riddle by identifying the location based on the given clues
    and prints the answer.
    """
    # Clue from the riddle
    distance_clue_km = 500

    # The identified answer and its details
    town_name = "Hanga Roa"
    island_name = "Easter Island (Rapa Nui)"
    bay_name = "Hanga Roa Bay"

    # Print the explanation and the final answer
    print(f"The island town more than {distance_clue_km} kilometers from another inhabited island is {town_name}.")
    print(f"{town_name} is the main town on {island_name}.")
    print(f"It sits on a bay called {bay_name}, which shares the town's name and is located on the volcanic island.")

solve_geography_riddle()