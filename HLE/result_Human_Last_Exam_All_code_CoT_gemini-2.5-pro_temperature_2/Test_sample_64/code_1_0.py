def find_pitch_accent():
    """
    This function identifies and explains the pitch accent for the Japanese word 「弟」.
    """
    word_kanji = "弟"
    word_reading = "おとうと (otouto)"
    correct_pattern_name = "Heiban"
    pitch_notation = "[0]"

    # Dictionary explaining the different pitch accent patterns
    patterns = {
        "Heiban": "The pitch starts low on the first mora, rises on the second, and remains high to the end of the word. A following particle is also high.",
        "Atamadaka": "The pitch is high only on the first mora and drops from the second mora onwards.",
        "Nakadaka": "The pitch rises after the first mora and falls on a specific mora in the middle of the word.",
        "Odaka": "The pitch rises after the first mora and stays high until the last mora. A following particle is pronounced low, causing a drop after the word."
    }

    # Print the result
    print(f"Word: {word_kanji} ({word_reading})")
    print("-" * 30)
    print(f"The standard pitch accent pattern is: {correct_pattern_name} ({pitch_notation})")
    print("\nExplanation:")
    print(f"The '{correct_pattern_name}' pattern means: {patterns[correct_pattern_name]}")
    
    # Illustrate the pattern for the specific word
    # おとうと has 4 morae: お(o), と(to), う(u), と(to)
    print("\nFor 「おとうと」, the pitch contour is:")
    print("Mora 1 (お): Low")
    print("Mora 2 (と): High")
    print("Mora 3 (う): High")
    print("Mora 4 (と): High")
    print("Pattern: L-H-H-H")

find_pitch_accent()