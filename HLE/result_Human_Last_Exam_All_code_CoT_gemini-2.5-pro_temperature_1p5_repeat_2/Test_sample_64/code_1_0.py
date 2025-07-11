import pyopenjtalk
import re

def get_pitch_accent_info(word):
    """
    Analyzes a Japanese word to determine its pitch accent pattern using pyopenjtalk.
    """
    try:
        # Get the detailed linguistic context label for the word.
        # It's a list, but for a single word, we only need the first element.
        full_context_label = pyopenjtalk.extract_fullcontext(word)[0]

        # Use regular expressions to find the accent position and mora count.
        # The accent position 'p3' is in the format /A:xx-xx+p3
        accent_match = re.search(r"/A:\d+-\d+\+(\d+)", full_context_label)
        # The mora count 'p5' is in the format /F:xx_p5
        mora_match = re.search(r"/F:\d+_(\d+)", full_context_label)

        if not accent_match or not mora_match:
            print(f"Could not parse accent information for '{word}'.")
            return

        accent_pos = int(accent_match.group(1))
        mora_count = int(mora_match.group(1))
        pronunciation = pyopenjtalk.g2p(word, kana=True)

        # Determine the accent pattern name based on the rules.
        accent_pattern = "Unknown"
        if accent_pos == 0:
            accent_pattern = "Heiban (A)"
        elif accent_pos == 1:
            accent_pattern = "Atamadaka (B)"
        elif accent_pos == mora_count:
            accent_pattern = "Odaka (D)"
        elif 1 < accent_pos < mora_count:
            accent_pattern = "Nakadaka (C)"

        # Print the detailed breakdown.
        print(f"Analysis for: 「{word}」")
        print(f"Pronunciation: {pronunciation}")
        print(f"Mora Count: {mora_count}")
        print(f"Accent Kernel Position: {accent_pos} (A value of 0 means there is no pitch drop within the word)")
        print(f"Resulting Pitch Accent Pattern: {accent_pattern}")

    except ImportError:
        print("Please install pyopenjtalk: pip install pyopenjtalk")
    except Exception as e:
        print(f"An error occurred: {e}")

# The word we want to analyze
target_word = "弟"
get_pitch_accent_info(target_word)