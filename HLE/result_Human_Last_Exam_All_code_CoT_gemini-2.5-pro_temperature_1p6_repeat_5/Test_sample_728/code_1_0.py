def solve_bechdel_test():
    """
    This function evaluates a predefined list of films against the Bechdel Test
    and prints the letters corresponding to the films that pass.
    """
    
    # Data structure representing the films and their Bechdel Test status.
    # The 'passes' boolean is True if a film meets all three criteria.
    # Data is based on public information from sources like bechdeltest.com.
    films_data = {
        'a': {'title': 'Igla', 'year': 1988, 'passes': True},
        'b': {'title': 'Ghost Town', 'year': 2008, 'passes': False},
        'c': {'title': 'Girls Will Be Girls', 'year': 2003, 'passes': True},
        'd': {'title': 'War Dogs', 'year': 2016, 'passes': False},
        'e': {'title': 'Slither', 'year': 2006, 'passes': True},
        'f': {'title': 'John Dies at the End', 'year': 2012, 'passes': True},
        'g': {'title': 'Man Who Knew Too Much, The', 'year': 1934, 'passes': False},
        'h': {'title': 'Ladies In Retirement', 'year': 1941, 'passes': True},
        'i': {'title': 'Internship, The', 'year': 2013, 'passes': False},
        'j': {'title': 'Tinker Bell and the Lost Treasure', 'year': 2009, 'passes': True}
    }

    # A list to hold the letters of films that pass the test.
    passing_films = []

    # Iterate through the films data and check the 'passes' status.
    for letter, data in films_data.items():
        if data['passes']:
            passing_films.append(letter)
            
    # Sort the list alphabetically to ensure consistent output order.
    passing_films.sort()

    # Join the letters with a comma and print the final result.
    print(','.join(passing_films))

solve_bechdel_test()