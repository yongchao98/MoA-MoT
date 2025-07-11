def find_words_with_non_first_syllable_accent():
    """
    This function identifies and lists Russian words from a specific text
    that have more than one syllable and do not have the accent on the first syllable.
    The analysis of syllable count and accent is pre-determined for accuracy.
    """
    # The list of words that fit the criteria, in order of appearance.
    # The analysis (syllable count and stress placement) has been done manually.
    # Шашлык (Шаш-лы́к) - 2 syllables, stress on 2nd.
    # горелым (го-ре́-лым) - 3 syllables, stress on 2nd.
    # прибежал (при-бе-жа́л) - 3 syllables, stress on 3rd.
    # сосед (со-се́д) - 2 syllables, stress on 2nd.
    # дошёл (до-шёл) - 2 syllables, stress on 2nd.
    # его (е-го́) - 2 syllables, stress on 2nd.
    qualifying_words = ["Шашлык", "горелым", "прибежал", "сосед", "дошёл", "его"]

    # Join the words into a single comma-separated string for output.
    result = ", ".join(qualifying_words)
    print(result)

find_words_with_non_first_syllable_accent()