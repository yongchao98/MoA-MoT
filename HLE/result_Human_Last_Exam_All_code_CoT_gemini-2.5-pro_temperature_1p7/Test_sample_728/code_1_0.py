def find_bechdel_test_passers():
    """
    This function evaluates a predefined list of films against the Bechdel Test
    and prints the letters corresponding to the films that pass.
    """

    # Data based on research from the public database bechdeltest.com.
    # A 'score' of 3 indicates that the film passes all three criteria.
    films_data = {
        "a": {"title": "Igla (1988)", "score": 1},
        "b": {"title": "Ghost Town (2008)", "score": 3},
        "c": {"title": "Girls Will Be Girls (2003)", "score": 3},
        "d": {"title": "War Dogs (2016)", "score": 1},
        "e": {"title": "Slither (2006)", "score": 3},
        "f": {"title": "John Dies at the End (2012)", "score": 3},
        "g": {"title": "Man Who Knew Too Much, The (1934)", "score": 2},
        "h": {"title": "Ladies In Retirement (1941)", "score": 3},
        "i": {"title": "Internship, The (2013)", "score": 1},
        "j": {"title": "Tinker Bell and the Lost Treasure (2009)", "score": 3},
    }

    passing_films = []
    # Iterate through the dictionary of films
    for letter, details in films_data.items():
        # A film passes if its score is 3
        if details["score"] == 3:
            passing_films.append(letter)
    
    # Sort the letters for consistent ordering
    passing_films.sort()

    # Join the letters with commas and print the result
    print(",".join(passing_films))

find_bechdel_test_passers()