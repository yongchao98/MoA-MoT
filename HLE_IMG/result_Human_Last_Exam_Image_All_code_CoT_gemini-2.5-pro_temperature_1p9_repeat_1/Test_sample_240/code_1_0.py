def solve_manuscript_query():
    """
    This function provides the structured answer to the four questions about the manuscript.
    1. Vowel-marking system name.
    2. Seven-century time span of origin.
    3. The material it is written on.
    4. The verse number of a specific location.
    """
    vowel_system = "babylonian"
    time_span = "09-15"
    material = "parchment"
    
    # Verse identification details
    book_abbreviation = "ezk"
    chapter_number = "023"
    verse_number = "030"
    verse_citation = f"{book_abbreviation}.{chapter_number}:{verse_number}"
    
    # Constructing the final answer string
    # The output format is: name yy-zz material book.ccc:vvv
    final_answer = f"{vowel_system} {time_span} {material} {verse_citation}"
    
    print(final_answer)

solve_manuscript_query()