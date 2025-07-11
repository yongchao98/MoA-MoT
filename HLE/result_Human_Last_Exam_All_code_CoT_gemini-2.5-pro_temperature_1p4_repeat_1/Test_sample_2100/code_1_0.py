def solve_riddle():
    """
    This function solves the riddle by analyzing the provided clues.
    """
    # Clue 1: The poetic description by Andrei Voznesensky.
    poet = "Andrei Voznesensky"
    description = "a joy-discovering sailor"
    poem_question = "Where is your grave, even a mound?"

    # Clue 2: The historical fact about the grave's discovery.
    discovery_period = "the late 1980s"

    print("Step 1: Analyzing the provided clues.")
    print(f" - The man was the subject of a poem by {poet}.")
    print(f" - The poet wondered about the location of the man's lost grave.")
    print(f" - The grave was finally discovered in {discovery_period} (specifically, in 1991).")
    print("-" * 20)

    print("Step 2: Connecting the clues to a person.")
    print(" - The poem in question is 'Saga' by Voznesensky.")
    print(" - The poem is about Vitus Bering, a famous Danish explorer in Russian naval service.")
    print(" - Vitus Bering died in 1741 on an island that was later named after him.")
    print(" - For over 250 years, the exact location of his grave was unknown.")
    print(" - In 1991, a Danish-Russian expedition found and identified his remains on Bering Island, fitting the timeline.")
    print("-" * 20)

    # Determine the final answer
    last_name = "Bering"

    print("Step 3: State the final answer.")
    print(f"The last name of the man is: {last_name}")

# Execute the solver function
solve_riddle()