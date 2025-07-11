def solve_riddle():
    """
    This function solves the historical riddle by analyzing the provided clues.
    """
    # Clue 1: The poetic description and the poet.
    poet = "Andrei Voznesensky"
    description = "a joy-discovering sailor"
    print(f"Analyzing Clue 1: The person was called '{description}' by {poet}.")
    print("This points towards a specific work by Voznesensky, likely romanticizing a Russian explorer.")
    print("-" * 30)

    # Clue 2: The lost grave.
    quote = "Where is your grave, even a mound?"
    print(f"Analyzing Clue 2: The quote '{quote}' indicates the person's burial site was lost.")
    print("-" * 30)

    # Clue 3: The discovery period.
    discovery_period = "late 1980s"
    print(f"Analyzing Clue 3: The answer to the grave's location emerged in the {discovery_period}.")
    print("-" * 30)

    # Synthesizing the clues.
    print("Connecting the dots:")
    print(f"Andrei Voznesensky wrote the libretto for the famous Russian rock opera 'Juno and Avos'.")
    print("The opera is about the Russian diplomat and explorer Nikolai Rezanov.")
    print("The opera, which premiered in 1981, contains the famous lament for Rezanov's lost grave.")
    print("Rezanov died in Krasnoyarsk in 1807, and his grave was later destroyed. Spurred by the opera's popularity and the era of Glasnost, his burial site was rediscovered and a monument was erected in the late 1980s.")
    print("-" * 30)

    # Final Answer
    last_name = "Rezanov"
    print("Therefore, the man's last name is:")
    print(last_name)

solve_riddle()