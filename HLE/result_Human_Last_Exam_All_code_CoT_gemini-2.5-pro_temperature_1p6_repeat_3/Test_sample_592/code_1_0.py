def find_destination_city():
    """
    This function stores movie plot information and retrieves the
    US city Tommaso moves to in the movie Ovosodo.
    """
    # Storing movie plot details in a dictionary.
    movie_database = {
        'Ovosodo': {
            'year': 1997,
            'character_moves': {
                'name': 'Tommaso',
                'destination_city': 'Philadelphia',
                'destination_country': 'United States'
            }
        }
    }

    # Retrieve the information from the data structure.
    movie_title = 'Ovosodo'
    character_name = movie_database[movie_title]['character_moves']['name']
    city = movie_database[movie_title]['character_moves']['destination_city']

    # Print the final answer.
    print(f"In the 1997 movie {movie_title}, the character {character_name} moves to the US city of {city}.")

find_destination_city()