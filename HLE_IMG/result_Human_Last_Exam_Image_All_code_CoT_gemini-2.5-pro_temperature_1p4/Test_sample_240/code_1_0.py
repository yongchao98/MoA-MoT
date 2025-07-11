def solve_manuscript_query():
    """
    Solves the four-part query about the manuscript image.
    """
    vowel_system = "babylonian"
    time_span = "07-13"
    material = "parchment"
    
    # Verse identification: Esther 1:2
    book = "est"
    chapter = 1
    verse = 2
    
    # Format the verse according to BHS notation Xbook.yyy:xxx
    verse_string = f"{book}.{chapter:03d}:{verse:03d}"
    
    # Combine all parts into the final answer string
    final_answer = f"{vowel_system} {time_span} {material} {verse_string}"
    
    print(final_answer)

solve_manuscript_query()