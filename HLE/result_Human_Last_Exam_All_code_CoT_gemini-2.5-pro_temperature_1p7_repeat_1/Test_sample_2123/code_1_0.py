def solve_poet_riddle():
    """
    Solves the riddle by finding an English poet whose surname matches
    a word used to describe the atmosphere of Vienna's boulevards.
    """

    # List of famous English poets' surnames
    poets = [
        "Shakespeare", "Wordsworth", "Byron", "Shelley", "Keats",
        "Tennyson", "Browning", "Eliot", "Auden", "Blake", "Burns", "Gay"
    ]

    # Words describing the vibrant, festive, and lively atmosphere of Vienna's boulevards
    vienna_descriptors = [
        "grand", "wide", "lively", "festive", "imperial", "gay"
    ]

    # Find the poet whose surname is also a descriptor
    found_surname = None
    for poet in poets:
        if poet.lower() in vienna_descriptors:
            found_surname = poet
            break

    if found_surname:
        print(f"The descriptive word associated with Viennese boulevards is '{found_surname.lower()}'.")
        print(f"This matches the surname of the English poet John {found_surname}.")
        print(f"The surname is: {found_surname}")
    else:
        print("Could not find a matching poet for the riddle.")

solve_poet_riddle()