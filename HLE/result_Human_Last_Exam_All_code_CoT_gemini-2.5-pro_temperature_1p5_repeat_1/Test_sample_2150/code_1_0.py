def solve_trilogy_riddle():
    """
    This function solves a three-part riddle to identify three cities.
    """

    # Clue 1: "In Irving Wallace's 'The Miracle', a young girl anticipates an
    # encounter with a lady in white. In what city does this scene take place?"
    # The novel is set in the French town famous for the Marian apparitions
    # witnessed by Bernadette Soubirous.
    city1 = "Lourdes"

    # Clue 2: "A poem by Yuri Krasnokutsky describes a cold stream flowing
    # through the blood of 'X'. In what city is 'X' located?"
    # The poem contains the line "Как будто бы в крови течёт холодная Нева"
    # which translates to "As if the cold Neva flows in the blood". The Neva river
    # flows through Saint Petersburg.
    city2 = "Saint Petersburg"

    # Clue 3: "The Croatian gastronomic paradise 'Dolac' has a nickname made
    # up of two words, one of which is a toponym. In what city was its
    # main analogue situated?"
    # Dolac market is in Zagreb and is nicknamed the "Belly of Zagreb". The most
    # famous "belly" of a city is "Les Halles" in Paris, nicknamed "le ventre de Paris"
    # (the belly of Paris) by Émile Zola.
    city3 = "Paris"

    # Print the final trilogy in the specified format.
    print(f"{{{city1}, {city2}, {city3}}}")

solve_trilogy_riddle()