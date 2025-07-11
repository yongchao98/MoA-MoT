def solve_art_history_query():
    """
    This function provides the answer to the user's question about the painting "The Radionist".
    """
    # Define the known variables based on the user's query and external knowledge.
    painting_title = "The Radionist"
    artist = "Kurt Günther"
    creation_year = 1927
    acquisition_year = 1967
    
    # The painting was acquired by the Galerie Neue Meister (New Masters Gallery),
    # which is part of the Staatliche Kunstsammlungen Dresden (Dresden State Art Collections).
    acquiring_museum = "Galerie Neue Meister, Staatliche Kunstsammlungen Dresden"

    # Print the final answer in a full sentence, including the numbers from the query.
    print(
        f'The {creation_year} tempera painting "{painting_title}" by {artist} '
        f'was acquired in {acquisition_year} by the {acquiring_museum}.'
    )

solve_art_history_query()