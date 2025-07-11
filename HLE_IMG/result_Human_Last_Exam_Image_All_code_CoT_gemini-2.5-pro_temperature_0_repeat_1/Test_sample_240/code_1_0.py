def generate_manuscript_analysis():
    """
    This function analyzes the provided manuscript fragment and generates the required information.
    """

    # 1. Vowel-marking system: The system is Babylonian, identified by the supralinear vowel points.
    vowel_system = "babylonian"

    # 2. Time span: The manuscript is dated to the 10th-12th c. A seven-century span
    #    from the 8th to the 14th century (08-14) correctly covers this period.
    start_century = 8
    end_century = 14
    time_span = f"{start_century:02d}-{end_century:02d}"

    # 3. Material: The manuscript is written on vellum, a type of parchment.
    material = "vellum"

    # 4. Verse Numbering: The first verse in the left column of the right page is 2 Samuel 20:11.
    #    The format is Xbook.yyy:xxx.
    book_prefix = "2"
    book_abbreviation = "sa"
    chapter_number = 20
    verse_number = 11
    verse_reference = f"{book_prefix}{book_abbreviation}.{chapter_number:03d}:{verse_number:03d}"

    # Combine the answers into the final string, all in lowercase.
    final_answer = f"{vowel_system} {time_span} {material} {verse_reference}"

    print(final_answer)

generate_manuscript_analysis()