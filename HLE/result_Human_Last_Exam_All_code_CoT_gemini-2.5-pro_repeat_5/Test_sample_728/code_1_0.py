def solve_bechdel_test():
    """
    Analyzes a predefined list of films against the Bechdel Test criteria
    and prints the letters of the films that pass.
    """
    # Film data with their Bechdel test evaluation.
    # True means the film passes all three criteria.
    # False means it fails at least one criterion.
    films = [
        {'letter': 'a', 'title': 'Igla (1988)', 'passes': False},
        {'letter': 'b', 'title': 'Ghost Town (2008)', 'passes': True},
        {'letter': 'c', 'title': 'Girls Will Be Girls (2003)', 'passes': True},
        {'letter': 'd', 'title': 'War Dogs (2016)', 'passes': False},
        {'letter': 'e', 'title': 'Slither (2006)', 'passes': True},
        {'letter': 'f', 'title': 'John Dies at the End (2012)', 'passes': False},
        {'letter': 'g', 'title': 'Man Who Knew Too Much, The (1934)', 'passes': False},
        {'letter': 'h', 'title': 'Ladies In Retirement (1941)', 'passes': True},
        {'letter': 'i', 'title': 'Internship, The (2013)', 'passes': False},
        {'letter': 'j', 'title': 'Tinker Bell and the Lost Treasure (2009)', 'passes': True}
    ]

    passing_films_letters = []
    for film in films:
        if film['passes']:
            passing_films_letters.append(film['letter'])

    # Joining the list of letters with a comma to match the required output format.
    result = ",".join(passing_films_letters)
    print(result)

solve_bechdel_test()