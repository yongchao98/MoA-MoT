def solve_riddle():
    """
    This function pieces together the clues from the riddle to identify the historical figure.
    """
    # Clue 1: The source of the quote points to a Russian context.
    clue_poet = "Andrei Voznesensky (a prominent Russian poet)"

    # Clue 2: The description of the man.
    clue_description = "'a joy-discovering sailor' (an explorer/navigator)"

    # Clue 3: The central mystery of the poem.
    clue_mystery = "A lost grave, whose location was unknown for centuries"

    # Clue 4: The time the mystery was solved.
    clue_discovery = "The grave was found in the late 1980s / early 1990s"

    # Synthesizing the clues points to the explorer whose remains were found in 1991 on Bering Island.
    solution = "Bering"

    print("Deducing the answer from the provided clues:")
    print("---------------------------------------------")
    # The final equation, showing how each clue contributes to the answer.
    print(f"The clue from '{clue_poet}' about")
    print(f"a '{clue_description}' with")
    print(f"a '{clue_mystery}' that was solved when")
    print(f"'{clue_discovery}'")
    print("=============================================")
    print(f"Logically leads to the last name: {solution}")

solve_riddle()