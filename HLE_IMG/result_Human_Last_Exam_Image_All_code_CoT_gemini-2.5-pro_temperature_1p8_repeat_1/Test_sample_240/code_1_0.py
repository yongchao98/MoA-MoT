def solve_manuscript_questions():
    """
    This function provides the answers to the four questions about the manuscript fragment.
    """
    # 1. Vowel-marking system: Babylonian (supralinear pointing)
    vowel_system = "babylonian"

    # 2. Time span: The Babylonian system was used from roughly the 7th to the 13th century.
    time_span = "07-13"

    # 3. Material: The fragment is written on treated animal skin.
    material = "parchment"

    # 4. Verse numbering: The first verse in the left column on the right page is 2 Samuel 21:1.
    #    Book: 2 Samuel -> 2sam
    #    Chapter: 21 -> 021
    #    Verse: 1 -> 001
    verse_number = "2sam.021:001"

    final_answer = f"{vowel_system} {time_span} {material} {verse_number}"
    print(final_answer)

solve_manuscript_questions()