def find_accompanying_book():
    """
    This function identifies and prints the name of the book of manners
    preserved in the same manuscript as the Middle English romance 'Sir Launfal'.
    
    The romance 'Sir Launfal' is found in the manuscript Cotton Caligula A.ii.
    A review of the contents of this manuscript reveals several romances and
    a specific text on etiquette.
    """
    
    manuscript = "Cotton Caligula A.ii"
    sir_launfal = "Sir Launfal"
    book_of_manners = "The Boke of Curtasye"

    print(f"The romance '{sir_launfal}' is preserved in the manuscript known as {manuscript}.")
    print(f"Alongside it, in the same manuscript, is the book of manners called '{book_of_manners}'.")

find_accompanying_book()