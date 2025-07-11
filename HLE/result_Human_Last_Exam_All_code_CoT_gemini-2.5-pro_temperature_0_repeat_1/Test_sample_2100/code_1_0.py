def find_historical_figure():
    """
    This function identifies a historical figure based on a set of clues.
    """
    # Clue 1: The poet who wrote about him.
    poet = "Andrei Voznesensky"

    # Clue 2: The description from the poem.
    description = "a joy-discovering sailor"

    # Clue 3: The historical context of his grave.
    grave_status = "Lost for many years, then rediscovered in the late 1980s."

    # Analysis: These clues point to the Russian explorer and diplomat Nikolai Rezanov.
    # Voznesensky's famous rock opera "Juno and Avos" is about Rezanov's journey.
    # Rezanov died in 1807, and his grave in Krasnoyarsk was lost until its rediscovery
    # in the late 1980s.
    last_name = "Rezanov"

    print(f"The man called '{description}' by {poet} is Nikolai {last_name}.")
    print(f"His grave was rediscovered in the { 'late 1980s' }.")
    print(f"The last name is: {last_name}")

find_historical_figure()