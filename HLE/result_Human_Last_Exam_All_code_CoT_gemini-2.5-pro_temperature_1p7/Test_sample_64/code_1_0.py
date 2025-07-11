def explain_pitch_accent():
    """
    Explains and visualizes the pitch accent pattern for the Japanese word 「弟」.
    """
    word_kanji = "弟"
    word_kana = "おとうと"
    morae = ["お", "と", "う", "と"]
    pattern_name = "Heiban (平板)"
    pattern_code = "[0]"

    print(f"The word is 「{word_kanji}」, which is read as 「{word_kana}」.")
    print(f"The standard pitch accent for this word is known as {pattern_name}, indicated by the notation {pattern_code}.")
    print("\nIn a Heiban (flat) pattern:")
    print("1. The first mora is pronounced with a low pitch.")
    print("2. The pitch rises on the second mora and stays high for the rest of the word.")
    print("3. A particle following the word will also maintain the high pitch.")

    # Generate visual representation of the pitch
    pitch = ["L"] + ["H"] * (len(morae) - 1)
    
    # Create the output string
    visual_pattern_word = []
    for i in range(len(morae)):
        visual_pattern_word.append(f"{morae[i]}({pitch[i]})")
        
    particle = "が"
    particle_pitch = "H"
    
    print("\nVisual Pitch Pattern:")
    print(f"Word by itself: {' '.join(visual_pattern_word)}")
    print(f"Word with a particle: {' '.join(visual_pattern_word)} {particle}({particle_pitch})")
    
    print("\nConclusion: The correct pattern is Heiban.")

explain_pitch_accent()