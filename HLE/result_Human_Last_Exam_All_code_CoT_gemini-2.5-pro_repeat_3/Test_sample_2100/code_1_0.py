def solve_riddle():
    """
    This function analyzes the clues to identify the historical figure.
    """
    # Clues from the prompt
    poet = "Andrei Voznesensky"
    description = "a joy-discovering sailor"
    grave_mystery = "Where is your grave, even a mound?"
    time_clue = "The answer emerged in the late 1980s"

    print("Analyzing the riddle's clues:")
    print(f"1. The figure was the subject of a poem by {poet}.")
    print(f"2. He was described as a '{description}'.")
    print(f"3. His grave was lost, prompting the question: '{grave_mystery}'?")
    print(f"4. The mystery surrounding his grave was revisited in the {time_clue}.")
    print("-" * 20)

    print("Connecting the clues:")
    print(f"The poet {poet} wrote the famous rock opera 'Juno and Avos'.")
    print("The main character of this opera is the Russian explorer and statesman Nikolai Rezanov.")
    print("Voznesensky's work portrays Rezanov as a visionary explorer and focuses on the tragedy of his death and lost grave in Siberia.")
    print("Interest in Rezanov's story and the search for his grave surged during the 'Glasnost' period of the late 1980s.")
    print("-" * 20)

    # The final answer
    last_name = "Rezanov"
    print(f"Based on the connection to Voznesensky's work and the historical context, the man's last name is:")
    print(last_name)

solve_riddle()