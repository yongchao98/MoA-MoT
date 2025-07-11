def solve_manuscript_query():
    """
    This function provides the answers to the four questions about the manuscript fragment.
    """
    
    # 1. Vowel-marking system name (first word, lowercase)
    # The vowel markings are above the letters, which is characteristic of the Babylonian system.
    vowel_system = "babylonian"
    
    # 2. Seven-century time span (format: xx-yy)
    # The Babylonian vocalization system was in use roughly from the 7th to 13th centuries.
    century_start = 7
    century_end = 13
    time_span = f"{century_start:02d}-{century_end:02d}"
    
    # 3. Material (lowercase)
    # The material is prepared animal skin, i.e., parchment.
    material = "parchment"
    
    # 4. Verse numbering (format: Xbook.yyy:xxx)
    # The text is identified as the beginning of the book of Ezekiel.
    # BHS (Biblia Hebraica Stuttgartensia) standard.
    book = "ezk"
    chapter = 1
    verse = 1
    verse_id = f"{book}.{chapter:03d}:{verse:03d}"
    
    # Combine answers into a single string
    final_answer = f"{vowel_system} {time_span} {material} {verse_id}"
    
    print(final_answer)

solve_manuscript_query()