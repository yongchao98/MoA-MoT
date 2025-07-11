def get_pitch_accent_info(word):
    """
    This function provides the standard pitch accent for a given Japanese word
    by looking it up in a predefined dictionary.
    """
    # Pitch accent database
    # [0] = Heiban, [1] = Atamadaka, etc.
    accent_data = {
        "弟": {
            "reading": "おとうと (otouto)",
            "pattern_type": "Heiban",
            "pattern_code": "[0]",
            "description": "The pitch starts low, rises on the second mora, and stays high.",
            "pitch_visualization": "お(↓) と(↑) う(→) と(→)"
        }
    }

    if word in accent_data:
        info = accent_data[word]
        print(f"Word: 「{word}」")
        print(f"Reading: {info['reading']}")
        print(f"Standard Pitch Accent Pattern: {info['pattern_type']} {info['pattern_code']}")
        print(f"Description: {info['description']}")
        print(f"Visualization: {info['pitch_visualization']}")
    else:
        print(f"Pitch accent information for 「{word}」 is not available.")

# The word in question
target_word = "弟"
get_pitch_accent_info(target_word)