def solve_riddle():
    """
    This function identifies and prints the answer to the geographical riddle.
    """
    # Key clues from the riddle
    clue_remoteness_km = 500
    clue_geography = "bay formed by a volcanic caldera"
    clue_naming = "town and bay share a name"

    # The identified location and its details
    town_name = "Hanga Roa"
    island_name = "Easter Island (Rapa Nui)"
    bay_name = "Hanga Roa Bay"
    distance_to_nearest_inhabited_land_km = 2075 # to Pitcairn Islands

    # Print the reasoning and the answer
    print("Searching for an island town with the following properties:")
    print(f"- More than {clue_remoteness_km} kilometers from another inhabited island.")
    print(f"- Sits on a {clue_geography}.")
    print(f"- The {clue_naming}.")
    print("\n--- Solution ---")
    print(f"The town is: {town_name}")
    print(f"It is located on: {island_name}")
    print(f"It sits on a bay named: {bay_name}")
    print(f"Fact check: The nearest inhabited land is over {distance_to_nearest_inhabited_land_km} km away.")
    print("Conclusion: Hanga Roa is the correct answer.")

solve_riddle()