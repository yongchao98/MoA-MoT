import unicodedata

def get_kanji_radical(character):
    """
    Looks up the radical for a given Kanji character.
    This is a simplified example with a predefined dictionary.
    """
    # A dictionary mapping Kanji to their radical's Hepburn name.
    # The key is the kanji and the value is the Hepburn transcription of the radical name.
    radical_map = {
        '謄': 'Gonben', # 言偏 (gonben)
        '媛': 'Onnahen', # 女偏 (onnahen)
        '松': 'Kihen',   # 木偏 (kihen)
        '海': 'Sanzui',  # さんずい (sanzui)
        '河': 'Sanzui',  # さんずい (sanzui)
        '聞': 'Mon',     # 門 (mon-gamae)
    }

    radical_name = radical_map.get(character)

    if radical_name:
        print(f"The radical for the character '{character}' is '{radical_name}'.")
        return radical_name
    else:
        print(f"Radical information for '{character}' not found in the dictionary.")
        return None

# The target character
target_character = '謄'
final_answer = get_kanji_radical(target_character)
# The final answer in the required format will be added after the code block.