def find_museum():
    """
    This function stores and prints the information about the acquisition
    of Kurt Günther's painting "The Radionist".
    """
    # Known information from the query
    artist = "Kurt Günther"
    painting_title = "The Radionist"
    creation_year = 1927
    acquisition_year = 1967
    
    # Information found through research
    acquiring_museum = "Lindenau-Museum Altenburg"
    
    # Print the final answer including all the numbers from the query
    print(f"The museum that acquired the {creation_year} tempera painting '{painting_title}' by {artist} in {acquisition_year} was the {acquiring_museum}.")

find_museum()