import fugashi
import sys

def get_pitch_accent(word):
    """
    This function looks up the pitch accent of a Japanese word using fugashi and unidic.
    """
    try:
        # Initialize the tagger. The unidic dictionary is required for accent info.
        tagger = fugashi.Tagger()
    except RuntimeError:
        print("Error: fugashi/MeCab could not be initialized.", file=sys.stderr)
        print("Please make sure you have installed fugashi and a dictionary:", file=sys.stderr)
        print("pip install fugashi unidic-lite", file=sys.stderr)
        sys.exit(1)

    # Process the word. We look at the first node as "弟" is a single word.
    node = tagger.parse(word)[0]

    # Extract relevant features from the dictionary entry.
    reading = node.feature.kana
    mora_count = len(reading)  # For オトート, length is 4, which is the mora count.
    
    # The accent type is provided as a string like "0", "1", etc.
    accent_type_str = node.feature.accent_type

    try:
        accent_pattern_num = int(accent_type_str)
    except (ValueError, TypeError):
        print(f"Could not determine accent pattern number for '{word}'.")
        return

    # Determine the pattern name based on the number and mora count.
    if accent_pattern_num == 0:
        accent_pattern_name = "Heiban"
        pitch_description = "L-H-H-H" # (L)ow-(H)igh
    elif accent_pattern_num == 1:
        accent_pattern_name = "Atamadaka"
        pitch_description = "H-L-L-L"
    elif 1 < accent_pattern_num < morae_count:
        accent_pattern_name = "Nakadaka"
        pitch_description = "L-H...H(fall)-L"
    elif accent_pattern_num == morae_count:
        accent_pattern_name = "Odaka"
        pitch_description = "L-H-H-H(fall after word)"
    else:
        accent_pattern_name = "Unknown"
        pitch_description = "N/A"

    print(f"Analysis for: 「{word}」")
    print(f"Reading (Kana): {reading}")
    print(f"Number of Morae: {mora_count}")
    print(f"Pitch Accent Number (from dictionary): {accent_pattern_num}")
    print(f"Pitch Pattern Type: {accent_pattern_name}")
    print(f"Pitch contour example for a {mora_count}-mora word of this type: {pitch_description}")
    
# The word we want to analyze.
target_word = "弟"
get_pitch_accent(target_word)
