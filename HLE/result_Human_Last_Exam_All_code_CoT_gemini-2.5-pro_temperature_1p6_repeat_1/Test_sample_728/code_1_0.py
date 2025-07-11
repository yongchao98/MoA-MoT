def solve_bechdel_test():
    """
    Evaluates a predefined list of films against the Bechdel Test and prints the
    letters corresponding to the films that pass.
    """
    # The Bechdel Test has three criteria:
    # 1. At least two named women.
    # 2. Who talk to each other.
    # 3. About something other than a man.

    # Data based on research from sources like bechdeltest.com.
    # The dictionary maps the film's letter to whether it passes all three criteria.
    film_data = {
        'a': {'title': 'Igla (1988)', 'passes': False},
        'b': {'title': 'Ghost Town (2008)', 'passes': False},
        'c': {'title': 'Girls Will Be Girls (2003)', 'passes': True},
        'd': {'title': 'War Dogs (2016)', 'passes': False},
        'e': {'title': 'Slither (2006)', 'passes': True},
        'f': {'title': 'John Dies at the End (2012)', 'passes': False},
        'g': {'title': 'Man Who Knew Too Much, The (1934)', 'passes': False},
        'h': {'title': 'Ladies In Retirement (1941)', 'passes': True},
        'i': {'title': 'Internship, The (2013)', 'passes': False},
        'j': {'title': 'Tinker Bell and the Lost Treasure (2009)', 'passes': True},
    }

    passing_films = []
    for letter, data in film_data.items():
        if data['passes']:
            passing_films.append(letter)

    # Sort the letters for a consistent output format
    passing_films.sort()

    # Print the final result as a comma-separated string
    print(','.join(passing_films))

solve_bechdel_test()