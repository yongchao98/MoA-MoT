def solve_bechdel_test():
    """
    This function evaluates a predefined list of films against the Bechdel Test
    and prints the letters of the films that pass.
    """
    # Film data based on research from bechdeltest.com and other movie databases.
    # The boolean value indicates if the film passes all three criteria.
    films = {
        'a': ("Igla", 1988, False),
        'b': ("Ghost Town", 2008, False),
        'c': ("Girls Will Be Girls", 2003, True),
        'd': ("War Dogs", 2016, False),
        'e': ("Slither", 2006, True),
        'f': ("John Dies at the End", 2012, False),
        'g': ("Man Who Knew Too Much, The", 1934, False),
        'h': ("Ladies In Retirement", 1941, True),
        'i': ("Internship, The", 2013, False),
        'j': ("Tinker Bell and the Lost Treasure", 2009, True)
    }

    passing_films = []
    for letter, data in films.items():
        if data[2]:  # Check if the film passes (the boolean is True)
            passing_films.append(letter)

    # Sort the letters alphabetically to ensure consistent output
    passing_films.sort()

    # Join the list of letters into a comma-separated string
    answer = ",".join(passing_films)
    print(answer)

solve_bechdel_test()