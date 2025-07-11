def solve_task():
    """
    This function provides the solution by printing the required information.
    1. The name of the author of the book.
    2. A list of the syllables with primary word stress for the first 10 words.
    """

    # Part 1: Author's name
    author_name = "maimonides"

    # Part 2: List of stressed syllables in Arabic script
    # The first 10 words are:
    # 1. أنا (ʾanā) -> stressed syllable is 'nā' (نا)
    # 2. أُبَيِّنُ (ʾubayyinu) -> stressed syllable is 'bay' (بي)
    # 3. لَكَ (laka) -> stressed syllable is 'la' (ل)
    # 4. أَيْضًا (ʾayḍan) -> stressed syllable is 'ʾay' (اي)
    # 5. فِي (fī) -> stressed syllable is 'fī' (في)
    # 6. هٰذَا (hādhā) -> stressed syllable is 'hā' (ها)
    # 7. الْفَصْلِ (al-faṣli) -> stressed syllable is 'faṣ' (فص)
    # 8. دَلَائِلَ (dalāʾila) -> stressed syllable is 'lā' (لا)
    # 9. التَّوْحِيدِ (al-tawḥīdi) -> stressed syllable is 'ḥī' (حي)
    # 10. عَلَىٰ (ʿalā) -> stressed syllable is 'lā' (لى)
    stressed_syllables = "نا بي ل اي في ها فص لا حي لى"

    print(author_name)
    print()  # For the single line space
    print(stressed_syllables)

solve_task()