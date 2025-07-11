def find_poet_and_artwork():
    """
    This function provides information about a specific Spanish poem.
    """
    poet = "Rafael Alberti"
    artwork = "the painting 'Guernica' by Pablo Picasso"
    collection = "'A la pintura' ('To Painting')"
    
    explanation = (
        f"The verses were written by the Spanish poet {poet}.\n"
        f"They were directed at {artwork}.\n\n"
        f"The lines are from a poem within Alberti's collection {collection}, "
        f"in which he describes the feeling of being trapped within the horrifying and chaotic scene depicted in the famous anti-war painting."
    )
    
    print(explanation)

find_poet_and_artwork()