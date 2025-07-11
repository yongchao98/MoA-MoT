def solve_bechdel_test():
    """
    Analyzes a predefined list of films against the Bechdel Test and
    prints the letters of the films that pass all three criteria.
    """
    # Film data and their Bechdel Test results based on research.
    # A film is marked 'True' if it passes all three conditions:
    # 1. It has at least two named women in it
    # 2. Who talk to each other
    # 3. About something besides a man
    film_data = {
        'a': {'title': 'Igla (1988)', 'passes': True},
        'b': {'title': 'Ghost Town (2008)', 'passes': False},
        'c': {'title': 'Girls Will Be Girls (2003)', 'passes': True},
        'd': {'title': 'War Dogs (2016)', 'passes': False},
        'e': {'title': 'Slither (2006)', 'passes': True},
        'f': {'title': 'John Dies at the End (2012)', 'passes': False},
        'g': {'title': 'Man Who Knew Too Much, The (1934)', 'passes': True},
        'h': {'title': 'Ladies In Retirement (1941)', 'passes': True},
        'i': {'title': 'Internship, The (2013)', 'passes': False},
        'j': {'title': 'Tinker Bell and the Lost Treasure (2009)', 'passes': True}
    }

    passing_films = []
    # Iterate through the film data to find which films pass
    for letter, data in film_data.items():
        if data['passes']:
            passing_films.append(letter)
            
    # Sort the list for consistent ordering (though it's already sorted)
    passing_films.sort()

    # Join the list into a comma-separated string for the final answer
    answer = ",".join(passing_films)
    
    print(answer)

solve_bechdel_test()