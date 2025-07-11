import collections

def find_pitch_accent():
    """
    Determines and explains the standard pitch accent for the Japanese word 「弟」.
    """
    word = "弟"
    reading = "おとうと"
    romaji = "otōto"
    morae = ["お", "と", "う", "と"]

    # In Japanese pitch accent notation, [0] indicates a Heiban (平板) pattern.
    # The standard accent for 「弟」(おとうと) is [0].
    accent_notation = 0
    correct_pattern_name = "Heiban"

    # Define the answer choices
    choices = {
        'A': 'Heiban',
        'B': 'Atamadaka',
        'C': 'Nakadaka',
        'D': 'Odaka',
        'E': 'Heiban or Nakadaka'
    }

    print(f"Analyzing the pitch accent for the word: 「{word}」 ({reading})")
    print("-" * 30)
    print(f"The word has {len(morae)} morae: {('-').join(morae)}")

    # Explanation of the Heiban pattern for this specific word
    print(f"The standard pitch accent for {reading} is the [{accent_notation}] pattern.")
    print(f"The [{accent_notation}] pattern is called 'Heiban' (平板).")
    print("In a Heiban pattern, the pitch starts low and rises on the second mora.")
    print("It then remains high through the end of the word and onto any following particle.")
    print("Pitch contour for おとうと: Low-High-High-High.")
    print("-" * 30)

    # Find the correct choice
    correct_choice = None
    for key, value in choices.items():
        if value == correct_pattern_name:
            correct_choice = key
            break

    print(f"Based on this analysis, the correct answer is A: {correct_pattern_name}.")

if __name__ == "__main__":
    find_pitch_accent()