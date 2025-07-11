def find_first_best_picture_with_luxor_obelisk():
    """
    This function provides the name of the first Academy Award Best Picture winner
    to feature a Luxor Obelisk on-screen.

    The logic is based on historical research of Best Picture winners and their filming locations.
    - The two Luxor Obelisks are in Luxor, Egypt, and Paris, France (Place de la Concorde).
    - Reviewing Best Picture winners chronologically, 'An American in Paris' is identified as the first.
    """
    
    # Information about the winning film
    movie_title = "An American in Paris"
    award_year = 1951
    obelisk_location_in_film = "Place de la Concorde, Paris"

    # Print the final answer detailing the film and the year.
    print(f"The first winner of the Academy Award for Best Picture to depict a Luxor Obelisk was:")
    print(f"Film: {movie_title}")
    print(f"Year of Award: {award_year}")
    print(f"The film shows the obelisk located in the {obelisk_location_in_film}.")

# Execute the function to print the result.
find_first_best_picture_with_luxor_obelisk()