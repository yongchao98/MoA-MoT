def solve_bechdel_test():
    """
    Identifies which films from a predefined list pass the Bechdel Test.
    """
    # The Bechdel test results for each film.
    # A film is marked True if it passes all three criteria.
    films = {
        'a': {"title": "Igla", "year": 1988, "passes": True},
        'b': {"title": "Ghost Town", "year": 2008, "passes": False},
        'c': {"title": "Girls Will Be Girls", "year": 2003, "passes": True},
        'd': {"title": "War Dogs", "year": 2016, "passes": False},
        'e': {"title": "Slither", "year": 2006, "passes": True},
        'f': {"title": "John Dies at the End", "year": 2012, "passes": True},
        'g': {"title": "Man Who Knew Too Much, The", "year": 1934, "passes": True},
        'h': {"title": "Ladies In Retirement", "year": 1941, "passes": True},
        'i': {"title": "Internship, The", "year": 2013, "passes": False},
        'j': {"title": "Tinker Bell and the Lost Treasure", "year": 2009, "passes": True},
    }

    passing_films_letters = []
    # Iterate through the dictionary and collect the letters of films that pass
    for letter, data in films.items():
        if data["passes"]:
            passing_films_letters.append(letter)
            
    # Sort the list alphabetically for consistent output
    passing_films_letters.sort()

    # Join the letters into a single comma-separated string and print
    print(','.join(passing_films_letters))

solve_bechdel_test()
<<<a,c,e,f,g,h,j>>>