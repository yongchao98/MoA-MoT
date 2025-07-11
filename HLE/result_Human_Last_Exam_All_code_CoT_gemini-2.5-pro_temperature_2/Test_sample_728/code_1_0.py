def solve_bechdel_test():
    """
    This function evaluates a predefined list of films against the Bechdel Test
    and prints the letters corresponding to the films that pass.
    """

    # Data based on research from bechdeltest.com
    # True = Passes all 3 criteria, False = Fails at least one criterion
    films = {
        'a': {'title': 'Igla (1988)', 'passes': True},
        'b': {'title': 'Ghost Town (2008)', 'passes': False},
        'c': {'title': 'Girls Will Be Girls (2003)', 'passes': True},
        'd': {'title': 'War Dogs (2016)', 'passes': False},
        'e': {'title': 'Slither (2006)', 'passes': True},
        'f': {'title': 'John Dies at the End (2012)', 'passes': False},
        'g': {'title': 'The Man Who Knew Too Much (1934)', 'passes': False},
        'h': {'title': 'Ladies In Retirement (1941)', 'passes': True},
        'i': {'title': 'The Internship (2013)', 'passes': False},
        'j': {'title': 'Tinker Bell and the Lost Treasure (2009)', 'passes': True}
    }

    passing_films = []
    for letter, data in films.items():
        if data['passes']:
            passing_films.append(letter)
    
    # Sort the letters alphabetically to ensure consistent output
    passing_films.sort()

    print(','.join(passing_films))

solve_bechdel_test()