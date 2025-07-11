def solve_riddle():
    """
    This function solves the geographical riddle based on its clues.
    """
    # The riddle provides a key number for the distance.
    distance_in_km = 500

    # Based on the clues, especially the extreme remoteness, we can identify the location.
    town_name = "Edinburgh of the Seven Seas"
    island_name = "Tristan da Cunha"
    geographical_feature = "a volcanic caldera"

    # Print the step-by-step reasoning, including the number from the riddle.
    print(f"The riddle describes an island town more than {distance_in_km} kilometers from another inhabited island.")
    print(f"This extreme isolation points to the island of {island_name}.")
    print(f"The town on this volcanic island is named '{town_name}'.")
    print(f"The island is a large volcano with {geographical_feature} at its summit, fitting the description.")
    print("\n---")
    print(f"Final Answer: {town_name}")

solve_riddle()