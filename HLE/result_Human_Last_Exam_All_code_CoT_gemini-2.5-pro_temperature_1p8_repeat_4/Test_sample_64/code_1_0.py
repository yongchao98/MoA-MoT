import collections

def get_pitch_accent(word):
    """
    Simulates looking up a word's pitch accent in a dictionary.
    Pitch accent patterns are denoted by a number:
    [0]: Heiban (平板) - Low-High-High...
    [1]: Atamadaka (頭高) - High-Low-Low...
    [2], [3], etc.: Nakadaka (中高) - Low-High...Low...
    [n] where n is number of mora: Odaka (尾高) - Low-High...High-Low(on particle)
    """
    pitch_accent_db = {
        '弟': {'kana': 'おとうと', 'pattern_num': 0, 'pattern_name': 'Heiban'},
        '頭': {'kana': 'あたま', 'pattern_num': 2, 'pattern_name': 'Nakadaka'},
        '猫': {'kana': 'ねこ', 'pattern_num': 1, 'pattern_name': 'Atamadaka'},
        '犬': {'kana': 'いぬ', 'pattern_num': 2, 'pattern_name': 'Nakadaka'},
        '花': {'kana': 'はな', 'pattern_num': 2, 'pattern_name': 'Odaka'},
    }

    if word in pitch_accent_db:
        info = pitch_accent_db[word]
        print(f"Word: 「{word}」")
        print(f"Kana: {info['kana']}")
        # We explicitly show the pattern number as part of the "equation" or result.
        print(f"The pitch accent pattern number is: {info['pattern_num']}")
        print(f"This corresponds to the pattern: {info['pattern_name']}")
        print("\nExplanation:")
        if info['pattern_name'] == 'Heiban':
            print("The pitch starts low, rises on the second mora, and then stays high through the rest of the word and onto any following particle (e.g., おとうとが -> L-H-H-H-H).")
    else:
        print(f"Pitch accent information for 「{word}」 not found in the database.")

# Look for the pitch accent of 「弟」
target_word = '弟'
get_pitch_accent(target_word)