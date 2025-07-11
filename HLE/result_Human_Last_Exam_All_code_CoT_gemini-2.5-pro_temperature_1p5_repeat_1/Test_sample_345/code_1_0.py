def solve_art_history_query():
    """
    This function provides the answer to the user's question about the acquisition
    of Kurt Günther's painting "The Radionist".
    """
    
    # Information from the query and external knowledge
    painting_year = 1927
    acquisition_year = 1967
    painting_title = "The Radionist"
    artist = "Kurt Günther"
    museum = "Lindenau-Museum Altenburg"

    # Constructing the sentence as an "equation" of facts, including each number.
    # Fact 1 (Painting Year): 1927
    # Fact 2 (Acquisition Year): 1967
    print(f"The painting '{painting_title}' from the year {painting_year} was acquired by the {museum} in the year {acquisition_year}.")

solve_art_history_query()