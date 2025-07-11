def solve_task():
    """
    This function provides the solution based on the analysis of the manuscript fragment.
    1. Vowel system: babylonian
    2. Time span: 07-13 (7th to 13th century)
    3. Material: parchment
    4. Verse reference: 2sam.004:004 (2 Samuel 4:4)
    """
    name = "babylonian"
    time_span = "07-13"
    material = "parchment"
    verse = "2sam.004:004"

    # The final output needs to be a single string with parts separated by spaces.
    # The verse numbering requires each number to be present.
    # For example, 2sam.004:004 needs to show all digits as requested.
    book_prefix = "2"
    book_name = "sam"
    chapter = "004"
    verse_num = "004"
    
    verse_string = f"{book_prefix}{book_name}.{chapter}:{verse_num}"

    final_answer = f"{name} {time_span} {material} {verse_string}"
    print(final_answer)

solve_task()