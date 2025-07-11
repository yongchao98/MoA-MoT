def get_movie_fact():
    """
    Addresses the user's question about the movie Ovosodo.
    """
    movie_title = "Ovosodo"
    character_name = "Tommaso"
    
    # Based on the film's plot, the character does not move to a US city.
    clarification = f"In the 1997 movie {movie_title}, the character {character_name} does not move to a city in the United States. The story is set in Livorno, Italy."
    
    print(clarification)

get_movie_fact()