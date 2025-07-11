def solve_trilogy():
    """
    This function solves the trilogy puzzle by deciphering each clue to find a city name.
    """

    # Clue 1: The novel "The Miracle" by Irving Wallace is about the apparitions
    # of Our Lady of Lourdes (the "lady in white") to a young girl,
    # which occurred in the French town of Lourdes.
    city1 = "Lourdes"

    # Clue 2: A well-known poem by Yuri Krasnokutsky contains the line
    # "холодный ручей струится в крови Москвы" which translates to
    # "a cold stream flows in the blood of Moscow". The city is Moscow.
    city2 = "Moscow"

    # Clue 3: The Dolac market in Zagreb is famously nicknamed the
    # "Belly of Zagreb" (Trbuh Zagreba). This is a direct reference to its
    # historical analogue, the Les Halles market in Paris, which was called
    # "The Belly of Paris" (Le Ventre de Paris).
    city3 = "Paris"

    # Print the final trilogy in the required format.
    print(f"{{{city1}, {city2}, {city3}}}")

solve_trilogy()