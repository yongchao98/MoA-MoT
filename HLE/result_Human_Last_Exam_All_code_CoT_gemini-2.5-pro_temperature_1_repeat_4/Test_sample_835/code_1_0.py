def solve_hebrew_syllables():
    """
    This function solves a two-part problem on Biblical Hebrew syllables and stress,
    and prints the formatted result.
    """

    # Part 1: Identify the primary stressed syllables in the first seven words of Psalm 74:1.
    # The primary stress is indicated by the main cantillation mark on the word.
    # The list contains the orthographic representation of each stressed syllable.
    stressed_syllables_q1 = [
        "כִּ֗יל",   # From מַשְׂכִּ֗יל (maśkîl), stress on the final syllable.
        "סָ֥ף",     # From לְאָ֫סָ֥ף (ləʾāsāp̄), stress on the final syllable marked by mercha.
        "מָ֣ה",     # From לָמָ֣ה (lāmāh), stress on the final syllable.
        "הִים",    # From אֱ֭לֹהִים (ʾĕlōhîm), stress is phonologically on the final syllable, despite the irregular placement of the tarcha accent.
        "נַ֣חְ",   # From זָנַ֣חְתָּ (zānaḥtā), stress on the penultimate syllable.
        "נֶ֑",      # From לָנֶ֑צַח (lāneṣaḥ), the word has penultimate stress on the syllable נֶ (ne). The etnachta (֑) accent is placed on the following syllable by rule, but marks the stress on נֶ.
        "שַׁ֥ן"      # From יֶעְשַׁ֥ן (yeʿšan), stress on the final syllable.
    ]
    answer1 = " ".join(stressed_syllables_q1)

    # Part 2: Identify the last syllable with marked secondary stress in 1 Chronicles 5:10.
    # Secondary stress is marked by a meteg ( ֽ ). The last occurrence in the verse is in the word לַגִּלְעָֽד׃ (laggilʿāḏ).
    # The meteg is on the final syllable, עָדֽ (ʿāḏ).
    answer2 = "עָדֽ"

    # Combine the answers according to the specified format: answer1,answer2
    final_answer = f"{answer1},{answer2}"

    print(final_answer)

solve_hebrew_syllables()