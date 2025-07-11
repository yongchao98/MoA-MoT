def find_pitch_accent():
    """
    Analyzes and prints the pitch accent pattern for the Japanese word 「弟」.
    """
    word_info = {
        "kanji": "弟",
        "reading": "おとうと",
        "pitch_pattern_number": 0,
        "morae_count": 4
    }

    kanji = word_info["kanji"]
    reading = word_info["reading"]
    pitch_number = word_info["pitch_pattern_number"]

    # In Japanese pitch accent notation:
    # [0] = Heiban (平板): Low-High-High-High...
    # [1] = Atamadaka (頭高): High-Low-Low-Low...
    # [2]...[N-1] = Nakadaka (中高): Pitch drops in the middle.
    # [N] = Odaka (尾高): Pitch is high until the end, but drops on a following particle.

    pattern_name = ""
    answer_choice = ""

    if pitch_number == 0:
        pattern_name = "Heiban (平板)"
        answer_choice = "A"
        explanation = (
            "A pitch number of [0] indicates a Heiban (flat) pattern. "
            "The pitch starts low on the first mora (お), rises on the second (と), "
            "and stays high for the rest of the word (うと)."
        )
    else:
        # This logic handles other cases, though not needed for this specific word.
        if pitch_number == 1:
            pattern_name = "Atamadaka (頭高)"
            answer_choice = "B"
        elif 1 < pitch_number < word_info["morae_count"]:
            pattern_name = "Nakadaka (中高)"
            answer_choice = "C"
        elif pitch_number == word_info["morae_count"]:
            pattern_name = "Odaka (尾高)"
            answer_choice = "D"
        explanation = "The word's pitch number does not match Heiban."


    print(f"Analysis for the word: 「{kanji}」 ({reading})")
    print(f"Standard pitch accent notation: [{pitch_number}]")
    print("-" * 30)
    print(f"Explanation: {explanation}")
    print(f"Conclusion: The pattern for 「{kanji}」 is {pattern_name}.")
    print(f"This corresponds to Answer Choice: {answer_choice}")


if __name__ == '__main__':
    find_pitch_accent()