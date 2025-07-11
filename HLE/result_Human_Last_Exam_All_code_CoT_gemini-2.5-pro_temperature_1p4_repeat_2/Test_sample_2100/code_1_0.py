def find_historical_figure():
    """
    Identifies a historical figure based on a set of clues.
    """
    # Clues provided in the query
    poet = "Andrei Voznesensky"
    description = "a joy-discovering sailor"
    grave_mystery_quote = "Where is your grave, even a mound?"
    time_frame = "late 1980s"

    # Analysis based on the clues
    # The poem is from the rock opera "Juno and Avos".
    # The opera is about the Russian diplomat and explorer Nikolai Rezanov.
    # Rezanov's grave in Krasnoyarsk was lost, and Voznesensky's work revived public interest in his fate,
    # effectively providing a cultural "answer" in the 1980s.
    last_name = "Rezanov"

    print(f"The man described by the poet {poet} is identified as Nikolai Rezanov.")
    print(f"The last name is: {last_name}")

find_historical_figure()