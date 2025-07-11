def find_bee_common_name_info():
    """
    This function provides the common name and year it was coined for the bee genus Perdita.
    The information is based on research of entomological literature.
    """
    # The currently used common name for the genus "Perdita" is "Fairy Bees".
    common_name = "Fairy Bees"

    # While the term may have been used colloquially earlier, a prominent, citable
    # use that helped popularize the name is in the 2014 book "The Bee: A Natural History".
    year_coined = 2014

    # The final output needs to be in the format "<name>, <year>".
    # The f-string below constructs this output.
    final_answer = f"{common_name}, {year_coined}"
    
    print(final_answer)

find_bee_common_name_info()