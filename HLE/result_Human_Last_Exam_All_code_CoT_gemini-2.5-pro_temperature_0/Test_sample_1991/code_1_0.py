def solve_task():
    """
    This function generates a fourteen-syllable phrase based on specific creative and formatting rules.
    """

    # The phrase is designed to have 14 syllables and describe abnormal cell death leading to cancer.
    # It is crafted to match the rhyme and polysyllabic nature of the 14-syllable Disney song title,
    # "Supercalifragilisticexpialidocious".
    #
    # Syllable breakdown:
    # A-po-pop-to-tic (5)
    # Dys-reg-u-la-tion (5)
    # Is (1)
    # A-tro-cious (3)
    # Total = 14 syllables.
    phrase = "Apopoptotic Dysregulation Is Atrocious"

    # The user requested that each word be capitalized.
    # Using .upper() ensures the entire phrase is capitalized for emphasis and to meet the requirement.
    final_phrase = phrase.upper()

    print(final_phrase)

solve_task()