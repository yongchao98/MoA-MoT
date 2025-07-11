def solve_bechdel_test():
    """
    This function evaluates a list of films against the Bechdel Test
    and identifies which ones pass all three criteria.
    """
    # Data representing the films and their Bechdel Test score (out of 3).
    # A score of 3 means the film passes all three rules.
    films = [
        {'letter': 'a', 'title': 'Igla', 'year': 1988, 'score': 1},
        {'letter': 'b', 'title': 'Ghost Town', 'year': 2008, 'score': 1},
        {'letter': 'c', 'title': 'Girls Will Be Girls', 'year': 2003, 'score': 3},
        {'letter': 'd', 'title': 'War Dogs', 'year': 2016, 'score': 1},
        {'letter': 'e', 'title': 'Slither', 'year': 2006, 'score': 3},
        {'letter': 'f', 'title': 'John Dies at the End', 'year': 2012, 'score': 1},
        {'letter': 'g', 'title': 'The Man Who Knew Too Much', 'year': 1934, 'score': 3},
        {'letter': 'h', 'title': 'Ladies In Retirement', 'year': 1941, 'score': 3},
        {'letter': 'i', 'title': 'The Internship', 'year': 2013, 'score': 1},
        {'letter': 'j', 'title': 'Tinker Bell and the Lost Treasure', 'year': 2009, 'score': 3},
    ]

    passing_films = []
    for film in films:
        # A film passes all three conditions if its score is 3.
        if film['score'] == 3:
            passing_films.append(film['letter'])

    # Format the output as a comma-separated string of lowercase letters.
    result = ",".join(passing_films)
    print(result)

solve_bechdel_test()
<<<c,e,g,h,j>>>