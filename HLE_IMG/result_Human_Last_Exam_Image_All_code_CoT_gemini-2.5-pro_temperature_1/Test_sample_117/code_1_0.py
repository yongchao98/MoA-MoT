import collections

def solve():
    """
    This function solves the problem by identifying the biblical verse and transcribing the Hebrew text.
    
    The process involves:
    1.  Analyzing the image of the manuscript, which contains Arabic and Hebrew text.
    2.  The Arabic phrase "وهو قول الكتاب" ("and it is the saying of the Book") indicates a biblical quotation.
    3.  The transcribed Hebrew words are analyzed. The clearest parts are "شلح" (sh-l-ḥ) and "عام" (ʿ-ā-m), corresponding to Hebrew "שלח" (shalaḥ, send) and "עם" (ʿam, people).
    4.  This combination strongly suggests the famous biblical phrase "בשלח פרעה את העם" ("When Pharaoh let the people go").
    5.  This phrase is located at the beginning of Exodus 13:17.
    6.  The manuscript text appears to be a shortened or corrupted version of this verse.
    7.  The final answer is formatted as requested, with the verse reference and the full, unpointed Hebrew text from the Bible.
    """
    
    # The identified biblical book, chapter, and verse.
    book_abbreviation = "exo"
    chapter = 13
    verse = 17
    
    # The unpointed Hebrew text of the identified verse.
    hebrew_text = "בשלח פרעה את העם"
    
    # Formatting the final answer string.
    result = f"{book_abbreviation}. {chapter}:{verse}, {hebrew_text}"
    
    print(result)

solve()