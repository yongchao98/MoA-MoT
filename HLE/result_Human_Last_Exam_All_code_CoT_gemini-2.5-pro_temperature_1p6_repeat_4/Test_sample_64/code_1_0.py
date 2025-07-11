import pandas as pd

def find_pitch_accent(word):
    """
    Simulates looking up a word's pitch accent in a dictionary.
    """
    # This data represents a small pitch accent dictionary.
    # The number [0] represents the Heiban pattern.
    accent_data = {
        '弟': {'reading': 'おとうと', 'accent_num': 0, 'pattern_name': 'Heiban'},
        '兄': {'reading': 'あに', 'accent_num': 1, 'pattern_name': 'Atamadaka'},
        '心': {'reading': 'こころ', 'accent_num': 2, 'pattern_name': 'Nakadaka'},
        '妹': {'reading': 'いもうと', 'accent_num': 4, 'pattern_name': 'Odaka'},
    }

    # Answer choices mapping
    choices = {
        'Heiban': 'A',
        'Atamadaka': 'B',
        'Nakadaka': 'C',
        'Odaka': 'D'
    }

    if word in accent_data:
        info = accent_data[word]
        print(f"Looking up the word: 「{word}」")
        print(f"Reading: {info['reading']}")
        print(f"Pitch Accent Number: {info['accent_num']}")
        print(f"This number corresponds to the '{info['pattern_name']}' pattern.")
        print(f"The correct answer choice is: {choices[info['pattern_name']]}")
    else:
        print(f"Pitch accent information for '{word}' not found.")

# Find the pitch accent for the target word.
find_pitch_accent('弟')