def perform_scansion():
    """
    This function contains the scansion for the given line of poetry.
    Line: "And in the letter, my cousin mentions a piece of advice"
    """
    
    # Syllable-by-syllable stress analysis:
    # And: Unstressed (x)
    # in: Unstressed (x)
    # the: Unstressed (x)
    # let-ter: Stressed, Unstressed (/) (x)
    # my: Unstressed (x)
    # cou-sin: Stressed, Unstressed (/) (x)
    # men-tions: Stressed, Unstressed (/) (x)
    # a: Unstressed (x)
    # piece: Stressed (/)
    # of: Unstressed (x)
    # ad-vice: Unstressed, Stressed (x) (/)
    
    scansion = "x" + "x" + "x" + "/" + "x" + "x" + "/" + "x" + "/" + "x" + "x" + "/" + "x" + "x" + "/"
    print(scansion)

perform_scansion()