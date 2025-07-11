def find_museum():
    """
    This function stores and prints the name of the museum that acquired
    the painting "The Radionist" in 1967. The information is based on
    data from art collection databases.
    """
    painting = "The Radionist"
    artist = "Kurt Günther"
    year_painted = 1927
    year_acquired = 1967
    
    # Based on provenance records, the acquiring museum is the Lindenau Museum Altenburg.
    museum = "Lindenau Museum Altenburg"

    print(f"The {year_painted} painting '{painting}' by {artist} was acquired in {year_acquired} by the following museum:")
    print(museum)

find_museum()