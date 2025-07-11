def find_accompanying_book():
    """
    This function provides the name of the book of manners
    preserved alongside the poem "Sir Launfal".
    """
    # The poem "Sir Launfal" is found in the manuscript Cotton Caligula A.ii.
    # This manuscript contains a collection of Middle English texts.
    manuscript = "Cotton Caligula A.ii"
    
    # Among the other works in this manuscript is a courtesy book, or a "book of manners".
    book_of_manners = "The Book of Curtesye"
    alternate_title = "Lytylle Childrenes Lytil Boke"

    print(f"The poem 'Sir Launfal' is preserved in the manuscript {manuscript}.")
    print(f"Alongside it is a book of manners known as '{book_of_manners}' (also called '{alternate_title}').")

find_accompanying_book()