import sys
import io

# Redirect stdout to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def find_pitch_accent(word):
    """
    Simulates looking up a word in a Japanese pitch accent dictionary.
    """
    # A small dictionary mapping words to their readings and accent patterns.
    # Pattern Notations:
    # [0] -> Heiban (平板): Starts low, rises, and stays high. L-H-H-H...
    # [1] -> Atamadaka (頭高): High on the first mora, then drops. H-L-L-L...
    # [n] -> Nakadaka (中高): Rise and fall are within the word. L-H-...-L...
    # [n] -> Odaka (尾高): Rises and stays high, but drops on the following particle. L-H-H... (L)
    pitch_accent_database = {
        "弟": {"reading": "おとうと", "pattern": "Heiban", "type": "A"},
        "雨": {"reading": "あめ", "pattern": "Atamadaka", "type": "B"},
        "心": {"reading": "こころ", "pattern": "Nakadaka", "type": "C"},
        "犬": {"reading": "いぬ", "pattern": "Odaka", "type": "D"}
    }

    if word in pitch_accent_database:
        info = pitch_accent_database[word]
        reading = info["reading"]
        pattern = info["pattern"]
        print(f"The word is 「{word}」, read as 「{reading}」.")
        print(f"The standard pitch accent pattern is: {pattern}")
    else:
        print(f"Information for the word '{word}' not found in the database.")

# The word we want to find the pitch accent for.
target_word = "弟"
find_pitch_accent(target_word)

# Restore stdout and print the captured output
sys.stdout = old_stdout
print(captured_output.getvalue())