def solve_hebrew_phonology():
    """
    This function solves the two-part Biblical Hebrew phonology task
    and prints the result in the specified format.
    """

    # Part 1: Analysis of Psalm 74:1 for primary stressed syllables.
    # The syllables are identified based on the placement of the main accent mark (ta'am).
    stressed_syllables_part1 = [
        "כִּ֗יל",  # from מַשְׂכִּ֗יל
        "אָ֫",     # from לְאָ֫סָ֥ף (stress marked by 'ole')
        "מָ֣ה",    # from לָמָ֣ה
        "לֹ֭",     # from אֱ֭לֹהִים
        "נַ֣ח",    # from זָנַ֣חְתָּ
        "נֶ֑",     # from לָנֶ֑צַח
        "שַׁ֥ן"     # from יֶעְשַׁ֥ן
    ]

    # Part 2: Analysis of 1 Chronicles 5:10 for the last secondary stressed syllable.
    # The syllable is identified by finding the last 'meteg' in the verse.
    # The last meteg is in the word לַגִּלְעָֽד׃ on the syllable עָֽ.
    secondary_stress_syllable_part2 = "עָֽ"

    # Combine the results as per the formatting instructions.
    # Part 1 syllables are joined by a space.
    part1_string = " ".join(stressed_syllables_part1)
    
    # Part 1 and Part 2 results are joined by a comma without a space.
    final_answer = f"{part1_string},{secondary_stress_syllable_part2}"

    print(final_answer)

solve_hebrew_phonology()