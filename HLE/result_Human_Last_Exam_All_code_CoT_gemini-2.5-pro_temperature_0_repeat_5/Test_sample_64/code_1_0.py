def find_pitch_accent(word):
    """
    This function simulates a lookup in a Japanese pitch accent dictionary.
    """
    accent_database = {
        "弟": {
            "reading": "おとうと",
            "morae": ["お", "と", "う", "と"],
            "pattern_name": "Heiban",
            "pitch_sequence": "Low-High-High-High",
            "answer_choice": "A"
        }
    }

    if word in accent_database:
        info = accent_database[word]
        print(f"Analysis for the word: 「{word}」")
        print(f"Reading: {info['reading']}")
        print(f"Morae: {' - '.join(info['morae'])}")
        print(f"Pitch Sequence: {info['pitch_sequence']}")
        print("-" * 30)
        print("Explanation:")
        print("The pitch starts low on the first mora 'お', rises on the second mora 'と', and remains high for the rest of the word.")
        print("This pattern is known as Heiban (平板).")
        print("-" * 30)
        print(f"Therefore, the correct answer choice is: {info['answer_choice']}")
    else:
        print(f"Pitch accent information for '{word}' not found.")

# The word we are looking for is 「弟」
target_word = "弟"
find_pitch_accent(target_word)