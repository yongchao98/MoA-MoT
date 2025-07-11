def find_unstressed_first_syllable_words():
    """
    This function identifies and lists Russian words from a specific text
    that have more than one syllable and do not have the accent on the first syllable.
    The analysis of syllables and accents is pre-determined for the words in the given sentence.
    """
    # The list of words that fit the criteria, based on linguistic analysis.
    # Words are in their order of appearance in the original text.
    # 1. Шашлык (шаш-лы́к) - 2 syllables, accent on the 2nd.
    # 2. горелым (го-ре́-лым) - 3 syllables, accent on the 2nd.
    # 3. прибежал (при-бе-жа́л) - 3 syllables, accent on the 3rd.
    # 4. сосед (со-се́д) - 2 syllables, accent on the 2nd.
    # 5. дошёл (до-шёл) - 2 syllables, accent on the 2nd.
    # 6. его (е-го́) - 2 syllables, accent on the 2nd.
    qualified_words = [
        "Шашлык",
        "горелым",
        "прибежал",
        "сосед",
        "дошёл",
        "его"
    ]

    # Join the words with a comma and print the result.
    result = ",".join(qualified_words)
    print(result)

find_unstressed_first_syllable_words()