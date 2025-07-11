def solve_riddle():
    """
    This function solves the geographical riddle by defining the criteria
    and printing the logical steps to reach the answer.
    """

    # --- Riddle Criteria ---
    # The prompt mentions a number that is part of the criteria.
    # Let's define it as a variable.
    min_distance_from_inhabited_land_km = 500

    # --- Logical Deduction ---
    print("Finding the island town based on the following criteria:")
    print(f"1. Must be more than {min_distance_from_inhabited_land_km} km from another inhabited island.")
    print("2. Sits on a bay formed by a volcanic caldera.")
    print("3. The town and the bay share the same name.")
    print("-" * 30)

    # Analysis of the clues
    analysis_step_1 = "The distance criterion points to a very remote location, like Easter Island (Rapa Nui)."
    analysis_step_2 = "Easter Island is volcanic, and its main town is Hanga Roa."
    analysis_step_3 = "In the Rapa Nui language, 'Hanga' means 'Bay' and 'Roa' means 'Long'."
    analysis_step_4 = "Therefore, the town's name, 'Hanga Roa', literally means 'Long Bay', satisfying the shared name criterion."

    # The final answer
    town_name = "Hanga Roa"
    island_name = "Easter Island (Rapa Nui)"

    print("Deduction:")
    print(f"- {analysis_step_1}")
    print(f"- {analysis_step_2}")
    print(f"- {analysis_step_3}")
    print(f"- {analysis_step_4}")
    print("-" * 30)
    print(f"The answer is the town of {town_name} on {island_name}.")

solve_riddle()