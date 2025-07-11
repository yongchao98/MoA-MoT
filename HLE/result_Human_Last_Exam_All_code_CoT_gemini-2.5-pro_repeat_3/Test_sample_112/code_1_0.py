def find_accompanying_book():
    """
    This function simulates a search through the contents of the manuscript
    Cotton MS Caligula A.ii to find the book of manners preserved alongside
    the romance 'Sir Launfal'.
    """
    
    # A representation of the major works in the manuscript
    manuscript_contents = {
        "Sir Launfal": "Romance",
        "The Boke of Curtasye": "Book of Manners",
        "Libeaus Desconus": "Romance",
        "Emaré": "Romance"
    }

    book_of_manners = None
    # Iterate through the dictionary to find the book of manners
    for title, genre in manuscript_contents.items():
        if genre == "Book of Manners":
            book_of_manners = title
            break

    if book_of_manners:
        print(f"The romance 'Sir Launfal' is preserved in a manuscript that also contains the following book of manners:")
        print(book_of_manners)
    else:
        print("Could not identify the book of manners in the manuscript's contents.")

find_accompanying_book()