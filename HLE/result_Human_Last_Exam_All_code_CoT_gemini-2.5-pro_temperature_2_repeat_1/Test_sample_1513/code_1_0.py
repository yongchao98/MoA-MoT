def solve_riddle():
    """
    Solves the geographical riddle by storing clues and presenting the answer.
    """
    # Clues from the riddle
    clue_remoteness_km = 500
    clue_location_type = "island town"
    clue_geography = "bay formed by a volcanic caldera"

    # Derived Answer based on analysis
    town_name = "Hanga Roa"
    island_name = "Easter Island (Rapa Nui)"
    bay_name = "Hanga Roa Bay"
    distance_to_nearest_inhabited_land_km = 2075 # To Pitcairn Islands

    # Explanation
    print(f"The riddle describes a very specific location.")
    print(f"Clue 1: An island town more than {clue_remoteness_km} km from another inhabited island.")
    print(f"Clue 2: It sits on a bay formed by a volcanic caldera sharing the town's name.")
    print("-" * 20)
    print(f"The location is {town_name} on {island_name}.")
    print(f"- It is extremely remote, over {distance_to_nearest_inhabited_land_km} km from the Pitcairn Islands.")
    print(f"- It sits on {bay_name}, on a volcanic island famous for its calderas.")
    print(f"Answer: {town_name}")

solve_riddle()