import sys

def find_kanji_radical(character):
    """
    Finds the Hepburn-romanized name of the radical for a given Kanji character.
    This function uses a predefined dictionary for the lookup.
    """
    # A dictionary mapping specific Kanji characters to their radical's name in Hepburn romanization.
    # The key is the kanji character, and the value is the name of its radical.
    radical_database = {
        '謄': 'Gonben',  # The target character. Its radical is 言 (gon), called Gonben as a left-side component.
        '媛': 'Onnahen', # The example character from the prompt. Its radical is 女 (onna), called Onnahen as a left-side component.
        '話': 'Gonben',  # Another example with the same radical.
        '林': 'Kihen',   # Radical is 木 (ki) on the left side.
        '海': 'Sanzui'   # Radical is 氵 (sanzui).
    }

    if character in radical_database:
        return radical_database[character]
    else:
        return f"Radical information for '{character}' is not available in the local database."

# The character we want to analyze
target_character = '謄'

# Find the radical name
radical_name = find_kanji_radical(target_character)

# Print the result
print(f"The bushu / radical in the Japanese character '{target_character}' is: {radical_name}")
