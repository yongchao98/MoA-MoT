def get_pitch_accent(word):
    """
    A simple function to look up the pitch accent of a Japanese word
    from a predefined dictionary.
    """
    # This dictionary simulates a pitch accent database.
    # The number [0] represents Heiban, [1] Atamadaka, [n] Nakadaka, etc.
    pitch_accent_db = {
        '弟': {
            'reading': 'おとうと (o-to-u-to)',
            'pattern_number': 0,
            'pattern_name': 'Heiban (平板)',
            'explanation': 'Starts low, rises, and stays high. The following particle is also high.'
        },
        '兄': {
            'reading': 'あに (a-ni)',
            'pattern_number': 1,
            'pattern_name': 'Atamadaka (頭高)',
            'explanation': 'Starts high and drops. The following particle is low.'
        },
        '心': {
            'reading': 'こころ (ko-ko-ro)',
            'pattern_number': 2,
            'pattern_name': 'Nakadaka (中高)',
            'explanation': 'Starts low, rises, then drops after the second mora. The following particle is low.'
        }
    }

    if word in pitch_accent_db:
        info = pitch_accent_db[word]
        print(f"Looking up the pitch accent for: 「{word}」")
        print(f"Reading: {info['reading']}")
        # As requested, outputting the number associated with the pattern
        print(f"Pattern Number: [{info['pattern_number']}]")
        print(f"Pattern Type: {info['pattern_name']}")
        print(f"Description: {info['explanation']}")
        print("\nBased on this information, the correct answer is Heiban.")
    else:
        print(f"Information for 「{word}」 not found.")

# The word we are interested in
target_word = '弟'
get_pitch_accent(target_word)