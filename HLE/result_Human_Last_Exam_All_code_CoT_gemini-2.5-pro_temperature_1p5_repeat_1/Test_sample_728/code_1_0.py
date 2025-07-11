def find_passing_bechdel_films():
    """
    Evaluates a predefined list of films against the Bechdel Test and
    prints the letters of the films that pass all three criteria.
    """

    # Data representing the films and their known Bechdel Test ratings.
    # The source for this data is bechdeltest.com.
    # A rating of 3 means the film passes all three criteria of the test.
    film_data = {
        'a': {'title': 'Igla', 'year': 1988, 'rating': 0},
        'b': {'title': 'Ghost Town', 'year': 2008, 'rating': 1},
        'c': {'title': 'Girls Will Be Girls', 'year': 2003, 'rating': 3},
        'd': {'title': 'War Dogs', 'year': 2016, 'rating': 0},
        'e': {'title': 'Slither', 'year': 2006, 'rating': 3},
        'f': {'title': 'John Dies at the End', 'year': 2012, 'rating': 1},
        'g': {'title': 'The Man Who Knew Too Much', 'year': 1934, 'rating': 3},
        'h': {'title': 'Ladies In Retirement', 'year': 1941, 'rating': 3},
        'i': {'title': 'The Internship', 'year': 2013, 'rating': 2},
        'j': {'title': 'Tinker Bell and the Lost Treasure', 'year': 2009, 'rating': 3}
    }

    # A list to store the letters of the films that pass the test
    passing_films_letters = []

    # Iterate over the film data
    for letter, data in film_data.items():
        # Check if the rating is 3, which signifies passing the test
        if data['rating'] == 3:
            passing_films_letters.append(letter)
    
    # Sort the letters alphabetically for a consistent output
    passing_films_letters.sort()

    # Join the list of letters into a single string, separated by commas
    result = ",".join(passing_films_letters)

    # Print the final result as requested
    print(result)

find_passing_bechdel_films()